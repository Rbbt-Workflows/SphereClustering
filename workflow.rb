
require 'rbbt-util'
require 'rbbt/workflow'

Workflow.require_workflow "COSMIC"
Workflow.require_workflow "Genomics"
Workflow.require_workflow "Structure"

module SphereClustering
  extend Workflow

  desc "Load mutations"
  input :protein, :string, "protein", ""
  input :mutations, :array, "Position	Samples", ""
  input :mutation_source, :select, "Select mutation source", "COSMIC", :select_options => ["mutations", "COSMIC"]
  task :load_mutations => :tsv do |protein, mutations, mutation_source|
    case mutation_source
    when "COSMIC"
      prot = protein.dup
      Protein.setup(prot, "Ensembl Protein ID", "Hsa")
      association_items = COSMIC.knowledge_base.subset(:mutation_protein_changes, :target => [prot], :source => :all)
      associations = association_items.tsv.to_double
      sample_mutations = COSMIC.knowledge_base.get_database(:sample_mutations, :type => :double, :merge => true, :target => "Sample name=~Sample", :source => "Genomic Mutation")
      sample_mutations.fields = ["Sample"]
      associations = associations.attach(sample_mutations)
      cosmic_mutations = TSV.setup({}, :key_field => "Position", :fields => ["Samples"], :type => :single)
      associations.each do |k,values|
        _genomic_mutation, _isoform, mutated_aa, samples = values


        mutated_aa = mutated_aa.first

        mutated_aa =~ /^([A-Z])(\d+)([A-Z])$/
        next if $1 == $3
        key = $2.to_i
        if cosmic_mutations.key? key
          cosmic_mutations[key] += samples.size
        else
          cosmic_mutations[key] = samples.size
        end
      end
      cosmic_mutations
    when "mutations"
      dumper = TSV::Dumper.new :key_field => "Position", :fields => ["Counts"], :type => :single, :into => dumper do |position,samples|
        [position,samples]
      end
      dumper
    end
  end


  input :protein, :string, "Ensembl Protein ID"
  input :distance, :integer, "Distance in Angstroms", 5
  input :pdbs, :array, "Only use these PDBs", []
  input :organism, :string, "Organism code", "Hsa/feb2014"
  task :protein_neighbourhood => :tsv do |protein,distance,pdbs,organism|

    protein = Protein.setup(protein.dup)
    ensp2seq = Organism.protein_sequence(organism).tsv :persist => true
    sequence = ensp2seq[protein]

    raise "#{protein.name} doesn't contains pdbs" if protein.pdbs.nil?
    pdbs_list = (pdbs.empty? ? protein.pdbs.keys : pdbs)

    uni = Organism.protein_identifiers(organism).index :target => "UniProt/SwissProt Accession", :persist => true
    uniprot = uni[protein]
    if uniprot and Structure::I3D_PROTEINS.include? uniprot
      filepos = Structure::I3D_PROTEINS.identify_field "FILENAME"
      Structure::I3D_PROTEINS[uniprot][filepos].each do |filename|
        type = filename =~ /EXP/ ? :pdb : :model
        url = "http://interactome3d.irbbarcelona.org/pdb.php?dataset=human&type1=proteins&type2=#{type}&pdb=#{filename}"
        pdb = filename.chomp(".pdb")
        if pdbs.empty? or pdbs_list.include?(pdb)
          pdbs_list.delete(pdb)
          pdbs_list.push(url)
        end
       end
    end

    tsv = TSV.setup({}, :key_field => "Sequence position", :fields => ["Neighbour residue"], :type => :flat, :cast => :to_i, :organism => organism)
    TSV.traverse pdbs_list, :bar => true, :type => :array, :cpus => 20, :into => tsv do |pdb|
      pdb_alignment = Structure.job(:pdb_alignment_map, [protein,pdb]*":", :sequence => sequence, :pdb => pdb).run
      next if pdb_alignment.size == 0

      pdb_neigbours = Structure.job(:neighbour_map, [protein, pdb]*":", :distance => distance, :pdb => pdb).run

      next if pdb_neigbours.size == 0

      pdb2seq = {}
      pdb_alignment.each do |k,values|
        values.each do |v|
          pdb2seq[v] = k
        end
      end

      result = []
      pdb_neigbours.each do |pos, ns|
        chain, _p = pos.split(":")
        res = pdb2seq[pos]
        next if res.nil?
        ns.each do |n|
          nchain, _np = n.split(":")
          next if nchain != chain
          nres = pdb2seq[n]
          next if nres.nil?

          result << [res, [nres]]
          result << [nres, [res]]
        end
      end
      result.extend MultipleResult

      result
    end

    tsv
  end

  input :protein, :string, "Ensembl Protein ID"
  input :spheres, :integer, "number of spheres", 3
  input :distance, :integer, "Distance in Angstroms", 5
  input :organism, :string, "Organism code", "Hsa/feb2014"
  dep :load_mutations, :protein => :protein
  dep :protein_neighbourhood, :protein => :protein, :distance => :distance
  task :protein_cluster => :tsv do |protein,spheres,distance,organism|
    prot = protein.dup
    Protein.setup(prot, "Ensembl Protein ID", "Hsa")

    neighbourhood = step(:protein_neighbourhood).load
    neighbourhood.each{|k,v| neighbourhood[k] = Set.new v}

    pos_freq_mutation_neighbourhood = {}

    pos_freq_mutation = step(:load_mutations).load

    total_cosmic_muts = pos_freq_mutation.map{|k,v| v.to_i}.reduce(:+)

    puts "Total cosmic mutations of the protein\t#{total_cosmic_muts}"
    pos_freq_mutation.each do |k,v|
      next if neighbourhood[k].nil?
      pos_freq_mutation_neighbourhood[k.to_i] = v.to_i
      neighbourhood[k].each do |n|
        next if pos_freq_mutation[n.to_s].nil?
        pos_freq_mutation_neighbourhood[k.to_i] += pos_freq_mutation[n.to_s].to_i
      end
    end

    sorted_positions = pos_freq_mutation_neighbourhood.sort_by{|p,c| c}.reverse.map{|p| [p[0].to_i, p[1]]}

    puts "Total of mutated positions of the protein\t#{sorted_positions.size}"
    positions = (0..spheres-1).to_a
    aa_size = sorted_positions.size-1
    raise "more spheres than mutated amino acids" if positions.size > aa_size
    candidates = []
    candidates.push([positions.map{|p| sorted_positions[p][1]}.reduce(:+), positions, 0])

    time = Time.now
    best = nil
    total_overlaps = 0
    overlaping_positions_hash = Hash.new(false)
    while not candidates.empty?
      best = candidates.pop

      _candidate_mutations, positions = best

      exist_overlap = false

      positions.size.times do |i|
        for j in i+1...positions.size
          key = "#{positions[i]}#{positions[j]}".to_i
          exist_overlap = overlaping_positions_hash[key]
          break if exist_overlap
          total_overlaps += 1
          first_pos = sorted_positions[positions[j]][0]
          second_pos = sorted_positions[positions[i]][0]
          exist_overlap = neighbourhood[first_pos.to_s].intersect? neighbourhood[second_pos.to_s]
          overlaping_positions_hash[key] = exist_overlap
          break if exist_overlap
        end
        break if exist_overlap
      end

      break unless exist_overlap
      best = nil

      for position in 0...positions.size
        condition = true
        position.times do |iteration|
          condition = condition && positions[iteration] == iteration
          break unless condition
        end
        next_to = ((position == positions.size-1) ? (aa_size) : (positions[position+1]))
        pos = ((next_to == aa_size) ? (positions[position]) : (positions[position]+1))
        condition = condition && (pos < next_to)
        if condition
          iterate_positions = positions.dup
          iterate_positions[position] += 1
          candidate = [iterate_positions.map{|p| sorted_positions[p][1]}.reduce(:+), iterate_positions]
          index = candidates.bsearch_index{|muts_positions| muts_positions[0] > candidate[0]}
          if index.nil?
            candidates << candidate
          else
            candidates.insert(index, candidate)
          end
        end
      end
    end

    if best.nil?
      "too much spheres"
    else

      c, position = best
      tsv = TSV.setup({}, :key_field => "Position", :fields => ["Mutations"], :type => :flat)
      positions.each_with_index do |p,i|
        pos, muts = sorted_positions[p]
        tsv[pos] ||= []
        tsv[pos] << muts
      end
      set_info :X, c

      tsv
    end
  end
  export_asynchronous :protein_cluster

  input :num_spheres, :array, "Sphere sizes", [1,2,3]
  input :distances, :array, "Distances", [5,7,9]
  dep :protein_cluster, :spheres => nil, :distance => nil, :compute => :bootstrap do |jobname,options|
    jobs = []
    options[:distances].each do |ns|
      options[:num_spheres].each do |d|
        jobs << SphereClustering.job(:protein_cluster, jobname, options.merge(:distance => d, :spheres => ns))
      end
    end
    jobs
  end

  task :calculate_pvalue => :tsv do
    table = TSV.setup({ })
    dependencies.each do |job|
      spheres = job.inputs[:spheres]
      distance = job.inputs[:distance]
      res =  job.info[:X]
    end
  end
end
