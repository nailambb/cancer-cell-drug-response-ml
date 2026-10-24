# read the CCLE expression matrix
load_expression(path)

# read GDSC and filter drug
load_response(path, drug_name)

# id mapping and tissue/lineage info
load_cell_line_metadata(path)

# map IDs and keep the intersection of cell lines
harmonize_ids(response_df, expression_df, metadata_df)

clean_features(X)

build_dataset(drug_name)

# write processed data into the data/processed folder
save_processed(...)


