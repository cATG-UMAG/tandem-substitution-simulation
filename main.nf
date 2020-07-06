#!/usr/bin/env nextflow

/*
 * Inputs preparation
 */
if(!(params.containsKey("groups")) || (params.groups == null)) {
  params.with_groups = false
} else {
  params.with_groups = true
}

references = file(params.references)
summary_config = file(params.summary_list)

// targets: [name, [...files]]
Channel.from(params.targets.entrySet())
  .map { [it.key, it.value.collect { file(it).baseName }] }
  .set { targets_grouped_ch }

// all targets combinations [file basename, name]
Channel.from(params.targets.entrySet())
  .flatMap { it.value.collect { x -> [file(x).baseName, it.key] } }
  .set { files_with_target_ch }

// all files from all targets in a flat list
Channel.from(params.targets.entrySet())
  .flatMap { it.value.collect { file(it) } }
  .into { targets_allfiles_ch; targets_allfiles_ch2 }

if (params.with_groups) {
  // groups: [name, [...targets]]
  Channel.from(params.groups.entrySet())
    .map { [it.key, it.value] }
    .set { groups_ch }

  // all groups combinations [target, name]
  Channel.from(params.groups.entrySet())
    .flatMap { it.value.collect { x -> [x, it.key] } }
    .into { targets_with_group_ch; targets_with_group_ch2 }
}


/*
 * Get observed tandems from fasta files
 */
process getObservedTandems {
  tag "$target"

  input:
  path fasta_file from targets_allfiles_ch
  path references

  output:
  tuple val(target), path("${target}.tsv") into observed_tandems_ch

  script:
  target = "${fasta_file.baseName}"
  """
  get_tandems.py $fasta_file $references ${target}.tsv
  """
}


/*
 * Join observed tandems under the same target
 */
targets_grouped_ch
 .flatMap { it[1].collect{ x -> [x, it[0]] } }
 .join(observed_tandems_ch)
 .collectFile(
      keepHeader: true,
      skip: 1,
      storeDir: 'output/tandem_info/observed/single/'
    ) {
    ["${it[1]}.tsv", it[2]]
 }
 .map { [it.baseName, it] }
 .set { observed_tandems_grouped_ch }


/*
 * Get mutational frequecies
 */
process getMutationalFrequencies {
  tag "$target"

  input:
  path sequences from targets_allfiles_ch2
  path references
  
  output:
  tuple val(target), path("${target}.tsv"), path("mbs_${target}.txt") into ungrouped_mutational_info_ch

  script:
  target = sequences.baseName
  """
  get_probabilities.py $sequences $references -m mbs_${target}.txt
  """
}


/*
 * Add target to mutational information and split channel
 * in the ones that need merging and the ones that not
 */
ungrouped_mutational_info_ch
  .join(files_with_target_ch)
  .groupTuple(by: 3)
  .map { [it[3], it[1], it[2]] }
  .branch {
    single: it[1].size() == 1
    multi: it[1].size() > 1
  }
  .set { multational_info_pre_ch }


/*
 * Merge mutational info when there are multiple files for a single target
 */
process mergeTargetMutationalInfo {
  tag "$name"
  publishDir 'output/mutation_info/single/', mode: 'copy', pattern: '*.tsv'

  input:
  tuple val(name), path(mutation_info_files), path(mutations_by_seq_files) from multational_info_pre_ch.multi

  output:
  tuple val(name), path("${name}.tsv"), path("mbs_${name}.txt") into merged_mutational_info_ch
  path "mbs_${name}.txt" into mutations_by_seq_ch

  script:
  """
  merge_mutation_info.py $mutation_info_files -o ${name}.tsv
  merge_mutations_by_seq.py $mutations_by_seq_files -n $name -o mbs_${name}.txt
  """
}


/*
 * Merge single and merged mutational info files in a single channel
 */
multational_info_pre_ch.single
  .map { it.flatten() }
  .concat(merged_mutational_info_ch)
  .multiMap {
    full1: it
    full2: it
    mutations_by_seq: it[2]
  }
  .set { mutational_info_ch }


/*
 * Merge all "mutations by sequence" files
 */
mutational_info_ch.mutations_by_seq
  .collectFile(
    name: 'mutations_by_seq.txt',
    storeDir: 'output/mutation_info/single/',
    newLine: true,
    sort: { it.text }
  )


/*
 * Make simulation
 */
process simulateTandems {
  tag "$name"
  cpus 4

  input:
  tuple val(name), path(mutation_info), path(mutations_by_seq) from mutational_info_ch.full1
  path references

  output:
  tuple val(name), path('tandems.tsv') into simulated_tandems_ch, simulated_tandems_ch2
  
  script:
  reference_names = params.targets[name].collect { file(it).baseName }.join(' ')
  productive_only = params.with_stopcodons ? "" : "-p"
  """
  tandem_generation.py $mutation_info $references $mutations_by_seq \
  ${params.n_simulations} tandems.tsv -r $reference_names -t ${task.cpus} $productive_only
  """
}


/*
 * Sort and compress tandem files
 */
process sortAndCompress {
  tag "$name"
  publishDir "output/tandem_info/simulated/single/", mode: 'copy'

  input:
  tuple val(name), path(tandems_file) from simulated_tandems_ch

  output:
  path "${name}.tsv.gz"

  script:
  """
  sort -nk1 -nk2 $tandems_file | gzip > ${name}.tsv.gz
  """
}


/*
 * Summarize tandem info
 */
if (params.with_sumaries) {
  process summarizeTandems {
    tag "$name"
    publishDir "output/tandem_info_summarized/single/", mode: 'copy'

    input:
    tuple val(name), path(tandems_file) from simulated_tandems_ch2
    path summary_config

    output:
    tuple val(name), path("${name}", type: 'dir') into summaries_ch

    script:
    """
    summarize_tandem_info.py $tandems_file $summary_config -o ${name}
    """
  }
}


/*
 * Grouping tasks: only executing them if there are groups defined
 */
if (params.with_groups) {
  /*
   * Merge tandem observed info by groups
   */
  groups_ch
    .flatMap { it[1].collect{ x -> [x, it[0]] } }
    .combine(observed_tandems_grouped_ch, by: 0)
    .collectFile(
        keepHeader: true,
        skip: 1,
        storeDir: 'output/tandem_info/observed/grouped/'
      ) {
      ["${it[1]}.tsv", it[2]]
  }

  /*
   * Arrange mutational info in groups to merge them later
   */
  mutational_info_ch.full2
    .combine(targets_with_group_ch, by: 0)
    .groupTuple(by: 3)
    .map { [it[3], it[1], it[2]] }
    .set { mutational_info_by_group_ch }


  /*
  * Merge mutational info by groups
  */
  process mergeGroupMutationalInfo {
    tag "$group_name"
    publishDir 'output/mutation_info/grouped/', mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(group_name), path(mutation_info_files), path(mutations_by_seq_files) from mutational_info_by_group_ch

    output:
    path "${group_name}.tsv"
    path 'mutations_by_seq.txt' into group_mutations_by_seq_ch

    script:
    """
    merge_mutation_info.py $mutation_info_files -o ${group_name}.tsv
    merge_mutations_by_seq.py $mutations_by_seq_files -n $group_name
    """
  }


  /*
  * Merge grouped "mutations by sequence" files
  */
  group_mutations_by_seq_ch
    .collectFile(
      name: 'mutations_by_seq.txt',
      storeDir: 'output/mutation_info/grouped/',
      newLine: true,
      sort: { it.text }
    )


  /*
  * Arrange summary directories in groups to merge them later
  */
  summaries_ch
    .combine(targets_with_group_ch2, by: 0)
    .groupTuple(by: 2)
    .map { [it[2], it[1]] }
    .set { summaries_by_group_ch }


  /*
  * Merge summaries by group
  */
  process mergeSummariesByGroup {
    tag "$group_name"
    publishDir "output/tandem_info_summarized/grouped/", mode: 'copy'

    input:
    tuple val(group_name), path(summary_dirs) from summaries_by_group_ch
    path summary_config

    output:
    path "${group_name}", type: 'dir'

    script:
    """
    merge_summaries.py $summary_config $summary_dirs -o "${group_name}"
    """
  }
}
