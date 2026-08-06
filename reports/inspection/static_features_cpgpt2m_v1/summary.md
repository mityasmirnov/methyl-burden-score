# CpGPT static feature export

- feature_set_id: `cpgpt2m_adapter_128_v1`
- genome_build: `GRCh38`
- source_commit: `a1d4f10d72cc30bdd6428b92ee6aa05e91adae21`
- checkpoint_sha256: `42b324bc3fd86062b00e2b3f742ae47c0012dabd2b480baaa3c6707ad12dc2f5`
- locus_table_sha256: `6b1358fd1a4f8c3b1ef0bf6d3faf7027c77d16e25dcd98d9f126053cba228be8`
- output_dimension: `128`
- storage_dtype: `float16`
- n_loci (registry): `1082522`
- n_mapped: `1076246`
- n_missing: `6276`
- mapping_rate: `0.994202`
- norm_mean: `10.131128`
- dim_var_mean: `0.838731`
- export_command: `uv run --extra cpgpt mbs features export-cpgpt --feature-set-id cpgpt2m_adapter_128_v1`
