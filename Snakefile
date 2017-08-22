configfile: "config.yaml"

rule all:
	input:
		expand('fastqc/{sample}_R1_fastqc.html', sample=config['samples']),
		expand('fastqc/{sample}_R2_fastqc.html', sample=config['samples']),
		expand('clean/{sample}_R1_paired.fastq.gz', sample=config['samples']),
		expand('clean/{sample}_R2_paired.fastq.gz', sample=config['samples']),
		expand('bam/{sample}.bam', sample=config['samples']),
		expand('bam/{sample}.bam.bai', sample=config['samples']),
		expand('bam/{sample}.bamqc', sample=config['samples']),
		expand('track/{sample}.tdf', sample=config['samples']),
		expand('peak_{p}/{treat}_peaks.bed', treat=config['treat'], p=config['p']),
		expand('peak_{p}/{treat}_peaks.xls', treat=config['treat'], p=config['p']),
		#expand('peak_{p}/{treat}_peaks_anno.xls', treat=config['treat'], p=config['p']),
		#expand('peak_{p}/{treat}_peaks_annoPie.pdf', treat=config['treat'], p=config['p']),
		#expand('figure/peak_{p}_correlation_heatmap_peak_score.pdf', p=config['p']),
		#expand('figure/peak_{p}_correlation_heatmap_read_count.pdf', p=config['p']),
		#expand('figure/peak_{p}_PCA.pdf', p=config['p']),
		#expand('RData/peak_{p}_dba.RData', p=config['p']),
		#expand('figure/peak_{p}_overlap.pdf', p=config['p']),
		## if has replicates, use diffbind 
		#expand('RData/peak_{p}_diff_test.RData', p=config['p']),
		#expand('RData/peak_{p}_diff_test_anno.RData', p=config['p']),
		#expand('table/peak_{p}_diff_test_anno.xlsx', p=config['p']),
		## if not, use MAnorm 
		#["MAnorm/{peak1}_vs_{peak2}_{p}/_all_peak_MAvalues.anno".format(peak1=x[0], peak2=x[1], p=config['p']) for x in zip(config['peak1'], config['peak2'])],
		#expand('table/all_peak_{p}_MAvalues_anno.xlsx', p=config['p']),
		#'doc/report.html'

rule fastqc:
	input:
		config['path']+'/{sample}_R1.fastq.gz',
		config['path']+'/{sample}_R2.fastq.gz'
	output:
		'fastqc/{sample}_R1_fastqc.html',
		'fastqc/{sample}_R2_fastqc.html'
	shell:
		'fastqc -t 2 -o fastqc {input}'

rule trimmomatic:
	input:
		r1 = config['path']+'/{sample}_R1.fastq.gz',
		r2 = config['path']+'/{sample}_R2.fastq.gz'
	output:
		r1_paired = 'clean/{sample}_R1_paired.fastq.gz',
		r2_paired = 'clean/{sample}_R2_paired.fastq.gz',
		r1_unpaired = 'clean/{sample}_R1_unpaired.fastq.gz',
		r2_unpaired = 'clean/{sample}_R2_unpaired.fastq.gz'
	params:
		adapter = config['adapter']
	shell:
		'trimmomatic PE -threads 3 -phred33 {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

rule bowtie2:
	input:
		r1 = 'clean/{sample}_R1_paired.fastq.gz',
		r2 = 'clean/{sample}_R2_paired.fastq.gz'
	output:
		bam = 'bam/{sample}.bam'
	params:
		prefix = 'bam/{sample}',
		cpu = config['cpu'],
		index = config['index']
	shell:
		"bowtie2 -p {params.cpu} -x {params.index} -1 {input.r1} -2 {input.r2}|samtools view -Sh -q 30 -F 4 -|grep -v 'XS:'|samtools view -Shub|samtools sort - -T {params.prefix} -o {output.bam}"

rule bam_idx:
	input:
		bam = 'bam/{sample}.bam'
	output:
		bai = 'bam/{sample}.bam.bai'
	shell:
		'samtools index {input.bam} {output.bai}'

rule bam_qc:
	input:
		bam = 'bam/{sample}.bam'
	output:
		bamqc = 'bam/{sample}.bamqc'
	params:
		cpu = config['cpu']
	shell:
		"qualimap bamqc --java-mem-size=10G -nt {params.cpu} -bam {input.bam} -outdir {output.bamqc}"

rule bam2bed:
	input:
		bam = 'bam/{sample}.bam'
	output:
		bed = 'MAnorm/{sample}.bed'
	shell:
		"bedtools bamtobed -split -i {input.bam} > {output.bed}"

rule bam2count:
	input:
		bam = 'bam/{sample}.bam'
	output:
		cnt = 'bam/{sample}.cnt'
	shell:
		"samtools view -c -F 4 {input.bam} > {output.cnt}"
		
rule bam2bedgraph:
	input:
		bam = 'bam/{sample}.bam',
		cnt = 'bam/{sample}.cnt'
	output:
		bg = 'track/{sample}.bedgraph'
	shell:
		"X=`awk '{{print 1/$1*1000000}}' {input.cnt}`; "
		"bedtools genomecov -ibam {input.bam} -bga -scale $X -split > {output.bg}"

rule bedgraph2tdf:
	input:
		bg = 'track/{sample}.bedgraph'
	output:
		tdf = 'track/{sample}.tdf'
	params:
		IGV = config['IGV']
	shell:
		"igvtools toTDF {input.bg} {output.tdf} {params.IGV}"

rule macs14:
	input:
		treat = 'bam/{treat}.chip.bam',
		control = 'bam/{treat}.input.bam'
	output:
		bed = 'peak_{p}/{treat}_peaks.bed',
		xls = 'peak_{p}/{treat}_peaks.xls'
	params:
		gsize = config['gsize'],
		pvalue = config['p'],
		prefix = 'peak_{p}/{treat}'
	shell:
		'macs14 -t {input.treat} -c {input.control} -g {params.gsize} -p {params.pvalue} -n {params.prefix} --nomodel'

rule peak_anno:
	input:
		bed = 'peak_{p}/{treat}_peaks.bed'
	output:
		xls = 'peak_{p}/{treat}_peaks_anno.xls',
		pie = 'peak_{p}/{treat}_peaks_annoPie.pdf'
	params:
		txdb = config['txdb'],
		gene2go = config['gene2go']
	shell:
		"Rscript script/peak_anno.R {input} {params} {output}"	

rule load_data_QC:
	input:
		['bam/{sample}.bam'.format(sample=x) for x in config['samples']],
		['peak_{p}/{treat}_peaks.bed'.format(p=config['p'], treat=x) for x in config['treat']]
	output:
		sheet = 'sample_sheet_peak_{p}.csv',
		heatmap_peak = 'figure/peak_{p}_correlation_heatmap_peak_score.pdf',
		heatmap_count = 'figure/peak_{p}_correlation_heatmap_read_count.pdf',
		PCA = 'figure/peak_{p}_PCA.pdf',
		RData = 'RData/peak_{p}_dba.RData'
	params:
		pvalue = config['p']	
	shell:
		"sed 's/$PEAK/peak_{params.pvalue}/' sample_sheet.csv > {output.sheet}; "
		"Rscript script/load_data_QC.R {output}"

rule peakset_overlap:
	input:
		'RData/peak_{p}_dba.RData'
	output:
		'figure/peak_{p}_overlap.pdf'
	params:
		mask = config['olp.mask']
	shell:
		'Rscript script/peakset_overlap.R {input} {params.mask} {output}'

rule diff_test:
	input:
		'RData/peak_{p}_dba.RData'
	output:
		'RData/peak_{p}_diff_test.RData'
	params:
		mask = config['diff.mask']
	shell:
		'Rscript script/diff_test.R {input} {params.mask} {output}'

rule diff_test_anno:
	input:
		'RData/peak_{p}_diff_test.RData'
	output:
		RData = 'RData/peak_{p}_diff_test_anno.RData',
		xlsx = 'table/peak_{p}_diff_test_anno.xlsx',
	params:
		txdb = config['txdb'],
		gene2go = config['gene2go']
	shell:
		'Rscript script/diff_test_anno.R {input} {params.txdb} {params.gene2go} {output.RData} {output.xlsx}'

rule MAnorm:
	input:
		peak1 = 'peak_{p}/{peak1}_peaks.bed',
		peak2 = 'peak_{p}/{peak2}_peaks.bed',
		read1 = 'MAnorm/{peak1}.bed',
		read2 = 'MAnorm/{peak2}.bed'
	output:
		dir = 'MAnorm/{peak1}_vs_{peak2}_{p}',
		xls = 'MAnorm/{peak1}_vs_{peak2}_{p}/_all_peak_MAvalues.xls'
	shell:
		"rm -r {output.dir}; MAnormFast --p1 {input.peak1} --r1 {input.read1} --p2 {input.peak2} --r2 {input.read2} -o {output.dir}"

rule MAnorm_anno:
	input:
		'MAnorm/{peak1}_vs_{peak2}_{p}/_all_peak_MAvalues.xls'
	output:
		'MAnorm/{peak1}_vs_{peak2}_{p}/_all_peak_MAvalues.anno',
		'MAnorm/{peak1}_vs_{peak2}_{p}/_all_peak_MAvalues.annoPie.pdf'
	params:
		txdb = config['txdb'],
		gene2go = config['gene2go']
	shell:
		"Rscript script/peak_anno.R {input} {params} {output}"	

rule merge_MAnorm_table:
	input:
		["MAnorm/{peak1}_vs_{peak2}_{p}/_all_peak_MAvalues.anno".format(peak1=x[0], peak2=x[1], p=config['p']) for x in zip(config['peak1'], config['peak2'])]
	output:
		'table/all_peak_{p}_MAvalues_anno.xlsx'
	shell:
		"Rscript script/merge_MAnorm_table.R config.yaml MAnorm {output}"	

rule report:
	input:
		rmd = 'doc/report.Rmd'
	output:
		html = 'doc/report.html'
	params:
		peak1 = config['samples'][0],
		peak2 = config['samples'][1]
	shell:
		"snakemake --config samples='{params.peak1}' treat='{params.peak1}' --dag | dot -Tsvg > doc/workflow.svg; "
		"Rscript script/make_report.R {input}"
		