#!bin/bash
python scripts/bpnet_extract_modisco_logos.py --ic_threshold 0.08 --modisco_model_file modisco/seq_width1000-lr0.001-lambda100-n_dil_layers9-conv_kernel_size7-tconv_kernel_size7-filters64/I_WT_D6CM_counts/modisco.h5 --output_h5 modisco/seq_width1000-lr0.001-lambda100-n_dil_layers9-conv_kernel_size7-tconv_kernel_size7-filters64/I_WT_D6CM_counts/modisco_logos.h5
python scripts/bpnet_extract_modisco_logos.py --ic_threshold 0.08 --modisco_model_file modisco/seq_width1000-lr0.001-lambda100-n_dil_layers9-conv_kernel_size7-tconv_kernel_size7-filters64/I_WT_S3MN_counts/modisco.h5 --output_h5 modisco/seq_width1000-lr0.001-lambda100-n_dil_layers9-conv_kernel_size7-tconv_kernel_size7-filters64/I_WT_S3MN_counts/modisco_logos.h5