import sys
import pysam

input_vcf = sys.argv[1]
output_vcf = sys.argv[2]

original_vcf = pysam.VariantFile(input_vcf)
formatted_vcf = pysam.VariantFile(output_vcf, "wz", header=original_vcf.header)

#assigned as tuple (Number=A)
for record in original_vcf:
    for sample in record.samples.values():
        if "VAF" not in sample:
            # VAF absent from FORMAT; force it to 0.0 rather than skipping
            try:
                sample["VAF"] = (0.0,)
            except Exception:
                pass  # if not in header, skip 
            continue
        vaf = sample["VAF"]
        if vaf is None:
            sample["VAF"] = (0.0,)
        #if vaf is . or non-empty
        else:
            sample["VAF"] = tuple(0.0 if x is None else x for x in vaf)
    formatted_vcf.write(record)

original_vcf.close()
formatted_vcf.close()