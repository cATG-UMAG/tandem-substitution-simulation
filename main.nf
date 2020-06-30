/*
 * Channels
 */
references = file(params.references)

Channel.fromPath(params.fasta_files)
  .into { fasta_files_ch; fasta_files_ch2 }

Channel.from((params.groups).entrySet())
  .map { [it.key, it.value] }
  .set { groups_ch }


/*
 * Get real tandems from fasta files
 */
process getRealTandems {
  tag "$target"
  publishDir "output/tandem_info/single/real/", mode: 'copy'

  input:
  path fasta_file from fasta_files_ch
  path references

  output:
  path "${target}.tsv" into real_tandems_ch

  script:
  target = "${fasta_file.baseName}"
  """
  get_tandems.py $fasta_file $references ${target}.tsv
  """
}


/*
 * Merge tandem real info into grouped files
 */
process mergeGroupTandemRealInfo {
  tag "$group_name"
  publishDir "output/tandem_info/grouped/real/", mode: 'copy'

  input:
  tuple val(group_name), val(group_items) from groups_ch
  path real_tandems_files from real_tandems_ch.collect()

  output:
  path "${group_name}.tsv"
  
  script:
  group_files = real_tandems_files.findAll { it.baseName in group_items }
  """
  head -n 1 ${group_files[0]} > ${group_name}.tsv &&
  for f in ${group_files.join(' ')}; do tail -q -n +2 \$f >> ${group_name}.tsv; done
  """
}


/*
 * Get mutational probabilities
 */
process getMutationProbabilities {
  tag "$target"
  publishDir "output/mutation_info/single/", mode: 'copy', pattern: "*.tsv"

  input:
  path sequences from fasta_files_ch2
  path references
  
  output:
  path "${target}.tsv" into mutational_frequencies_ch
  path "mutations_by_seq.txt" into mutations_by_seq

  script:
  target = sequences.baseName
  """
  get_probabilities.py $sequences $references
  """
}


// merge mutations_by_seq files
mutations_by_seq
  .collectFile(
    name: "mutations_by_seq.txt",
    storeDir: "output/mutation_info/single/",
    newLine: true,
    sort: { it.text }
  )