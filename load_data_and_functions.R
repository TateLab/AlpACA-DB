## Explicitly install required packages

library("tidyverse")
library("r3dmol")
library("fs")
library("data.table")
library("shinythemes")
library("bio3d")
library("shiny")
library("splitstackshape")
library("RSQLite")
library("DBI")
library("dbplyr")
library("dplyr") ## To allow for SQLite queries -- issues if loaded before




### Use SQLite database of res.af -- minimise RAM for Shiny app online
sqlite<-dbDriver("SQLite")

## pPSE SQL
pPSE_connect <- dbConnect(sqlite,"shiny_data/pPSE.sqlite")
res.af <- dplyr::tbl(pPSE_connect, "pPSE_df")

## Ligandability table
ligandability_connect <- dbConnect(sqlite,"shiny_data/ligandability_db.sqlite")
ligandability_df <- dplyr::tbl(ligandability_connect, "ligandability_df")


## Uniprot ID to gene name matching
if(sum(str_detect(ls(), "uniprot_id_gene_match")) == 0){
  ### Table with uniprot-to-gene name mapping
  uniprot_id_gene_match <- dir_ls("shiny_data/resources", regexp = "FASTA") %>%
    fread() %>%
    mutate(protein_id = sub(Gene, 
                            pattern = "..\\|(.*)\\|.*",
                            replacement = "\\1")) %>%
    select(protein_id, Genes) %>%
    rename(gene_name = "Genes")
  
}


## Return gene name from uniprot ID
return_gene_name <- function(input_uniprot_id){
  
  if(input_uniprot_id != ""){
    
    temp.gene_name <- uniprot_id_gene_match %>%
      filter(protein_id == input_uniprot_id) %>%
      pull(gene_name)
    
    paste0("AlphaFold structure for ", 
           return_gene_name(temp.gene_name))
    
  }
  
  return(temp.gene_name)
  
}

## Divide each function up - with final function combining all


# 1. Download AlphaFold structure
download_alphafold_prediction = function(input_uniprot_id, 
                                         input_dir){
  
  ### Check if directory has '/' at end
  if(!str_detect(input_dir, ".*\\/$")){
    
    dir_to_download_to <- paste0(input_dir, "/", input_uniprot_id)
    
  }else{
    
    dir_to_download_to <- paste0(input_dir, input_uniprot_id)
    
  }
  
  ### Create temp directory
  
  ## Only create new directory if not present
  if(!dir.exists(dir_to_download_to)){
    
    dir.create(dir_to_download_to)
    
  }
  
  # Download to folder
  utils::download.file(paste0("https://alphafold.ebi.ac.uk/files/AF-", 
                              input_uniprot_id,
                              "-F1-model_v4.pdb"), 
                       destfile = paste0(dir_to_download_to, "/AlphaFold_", input_uniprot_id, "_Original.pdb"))
  
  
}

# 2. Get path to downloaded .pdb structure
alphafold_download_dir <- function(input_uniprot_id, 
                                   input_dir){
  
  ### Check if directory has '/' at end
  if(!str_detect(input_dir, ".*\\/$")){
    
    dir_to_download_to <- paste0(input_dir, "/", input_uniprot_id)
    
  }else{
    
    dir_to_download_to <- paste0(input_dir, input_uniprot_id)
    
  }
  
  return(dir_to_download_to)
  
}


# 3. Reformat with pPSE values and visualise
pPSE_alphafold_visualise <- function(input_uniprot_id, 
                                     input_dir,
                                     start_col = "blue",
                                     end_col = "red"){
  
  dir_to_download_to <- alphafold_download_dir(input_uniprot_id,
                                               input_dir)
  
  ## Downloaded AlphaFold structure = original_af
  download_alphafold_prediction(input_uniprot_id,
                                input_dir)
  
  
  ## Re-formatted AlphaFold structure = reformat_af
  original_af <- dir_ls(dir_to_download_to, regexp = "AlphaFold.*_Original.pdb") %>%
    read.pdb()
  
  # Make table for specified uniprot_id from res.af
  temp.af.tbl <- res.af %>%
    filter(protein_id == input_uniprot_id) %>%
    as.data.frame()
  
  ## Get lengths of residue table from downlaoded AF structure
  temp.af.pdb.length <- original_af$atom %>%
    pull(resno) %>%
    unique() %>%
    length()
  
  ## Get length of pPSE table from res.af
  temp.af.tbl.length <- temp.af.tbl %>%
    nrow()
  
  
  if(temp.af.pdb.length == temp.af.tbl.length){  ## Only match if the lengths are the same
    
    ## Summarise number of rows per AA in AF structure (depends on individual AA #CNSO atoms)
    temp.af.natom.per.res <- original_af$atom %>%
      group_by(resno) %>%
      summarise(N.Rows=n()) %>%
      rename(position = "resno")
    
    temp.af.nres.join <- temp.af.tbl %>%
      as.data.frame() %>%
      left_join(as.data.frame(temp.af.natom.per.res)) %>%
      expandRows(count = "N.Rows", ## Replicate rows as specificed in AF structure ($atom table)
                 count.is.col = T)
    
    original_af$atom$b <- temp.af.nres.join$pPSE ## Replace B-factors with pPSE values
    
    modified_pdb_file <- paste0(dir_to_download_to, "/AlphaFold_", input_uniprot_id, "_Modified_pPSE.pdb")
    
    ## Write modified structure to temporary folder
    write.pdb(original_af,
              file = modified_pdb_file)
    
    max.pPSE <- temp.af.tbl %>%
      pull(pPSE) %>%
      max()
    
    ### Set colour gradient
    pPSE_gradient <- paste0('"', 
                            paste(colorRampPalette(c(start_col, end_col))(max.pPSE), 
                                  collapse = '", "'), 
                            '"')
    
    ### Plot re-formatted 3D structure
    r3dmol() %>%
      m_add_model(data = modified_pdb_file, format = "pdb") %>%
      m_set_style(style = m_style_cartoon(
        colorfunc = paste0("
              function(atom) {
                const color = [", pPSE_gradient,"]
                return color[Math.round(atom.b)]
              }")
      )) %>%
      m_zoom_to() 
    
  }
  
  
}


# 4. Use existing quality values
quality_alphafold_visualise <- function(input_uniprot_id, 
                                        input_dir,
                                        start_col = "blue",
                                        end_col = "red"){
  
  temp_dir <- alphafold_download_dir(input_uniprot_id,
                                     input_dir)
  
  dir_to_download_to <- paste0(temp_dir, "/AlphaFold_", input_uniprot_id, "_Original.pdb")
  
  
  ### Set colour gradient
  quality_gradient <- paste0('"', 
                             paste(colorRampPalette(c(start_col, end_col))(101), 
                                   collapse = '", "'), 
                             '"')
  
  ### Path for downlaoded AF structure
  
  ### Plot re-formatted 3D structure
  r3dmol() %>%
    m_add_model(data = dir_to_download_to, format = "pdb") %>%
    m_set_style(style = m_style_cartoon(
      colorfunc = paste0("
            function(atom) {
              const color = [", quality_gradient,"]
              return color[Math.round(atom.b)]
            }")
    )) %>%
    m_zoom_to() 
  
  
}


# 5. Colour specific residue(s) differently
specific_residue_visualise <- function(input_uniprot_id, 
                                       input_dir,
                                       residues, ## Accepts comma separated values or single set of continuous e.g. 1-10
                                       residue_col = "blue",
                                       rest_of_structure_col = "white"){
  ### Handle X numbers in format: 1,2,3,4 or 1-4
  if(str_detect(residues, "-")){
    
    get_start <- sub(residues, pattern = "(\\d*)-\\d*", replacement = "\\1") %>%
      as.numeric()
    get_end <- sub(residues, pattern = "\\d*-(\\d*)", replacement = "\\1") %>%
      as.numeric()
    
    reformat_residues <- get_start:get_end
    
  }else if(str_detect(residues, ",")){
    
    reformat_residues <- gsub(residues, pattern=" ", replacement="") %>%
      str_split(pattern=",", simplify = T) %>%
      as.numeric() %>%
      sort()
    
    
  }else{
    
    reformat_residues <- as.numeric(residues)
    
  }
  
  temp_residues <- reformat_residues
  
  dir_to_download_to <- alphafold_download_dir(input_uniprot_id,
                                               input_dir)
  
  ## Downloaded AlphaFold structure = original_af
  ## Re-formatted AlphaFold structure = reformat_af
  original_af <- dir_ls(dir_to_download_to, regexp = "AlphaFold.*_Original.pdb") %>%
    read.pdb()
  
  # Make temporary table to modify b values and put back into $atom table
  residue_edit_df <- original_af$atom %>%
    mutate(b = ifelse(resno %in% temp_residues, ## Accepts multiple residues
                      yes = 1, ## Arbitrary values to 
                      no = 0)) ## distinguish colours for plotting
  
  ## Replace b-factor column with new column based on specified residues
  original_af$atom$b <- residue_edit_df$b ## Replace B-factors with pPSE values
  
  ## Write out filename for modified pdb structure
  modified_structure_file <- paste0(dir_to_download_to, "/AlphaFold_", input_uniprot_id, "_Modified_SpecificResidue.pdb")
  
  ## Write modified structure to temporary folder
  write.pdb(original_af,
            file = modified_structure_file)  
  
  ### Set two colour options (specified residues, rest of structure)
  specific_residue_gradient <- paste0('"', 
                                      paste(colorRampPalette(c(rest_of_structure_col, 
                                                               residue_col))(2), 
                                            collapse = '", "'), 
                                      '"')
  
  ### Plot re-formatted 3D structure
  r3dmol() %>%
    m_add_model(data = modified_structure_file, format = "pdb") %>%
    m_set_style(style = m_style_cartoon(
      colorfunc = paste0("
              function(atom) {
                const color = [", specific_residue_gradient,"]
                return color[Math.round(atom.b)]
              }")
    )) %>%
    m_zoom_to() 
  
}


# 6. Highlight all Cys
all_cys_visualise <- function(input_uniprot_id, 
                              input_dir,
                              residue_col = "blue",
                              rest_of_structure_col = "white"){
  
  dir_to_download_to <- alphafold_download_dir(input_uniprot_id,
                                               input_dir)
  
  ## Downloaded AlphaFold structure = original_af
  ## Re-formatted AlphaFold structure = reformat_af
  original_af <- dir_ls(dir_to_download_to, regexp = "AlphaFold.*_Original.pdb") %>%
    read.pdb()
  
  # Make temporary table to modify b values and put back into $atom table
  residue_edit_df <- original_af$atom %>%
    mutate(b = ifelse(resid == "CYS", ## Accepts multiple residues
                      yes = 1, ## Arbitrary values to 
                      no = 0)) ## distinguish colours for plotting
  
  ## Replace b-factor column with new column based on specified residues
  original_af$atom$b <- residue_edit_df$b ## Replace B-factors with pPSE values
  
  ## Write out filename for modified pdb structure
  modified_structure_file <- paste0(dir_to_download_to, "/AlphaFold_", input_uniprot_id, "_Modified_AllCys.pdb")
  
  ## Write modified structure to temporary folder
  write.pdb(original_af,
            file = modified_structure_file)  
  
  ### Set two colour options (specified residues, rest of structure)
  specific_residue_gradient <- paste0('"', 
                                      paste(colorRampPalette(c(rest_of_structure_col, 
                                                               residue_col))(2), 
                                            collapse = '", "'), 
                                      '"')
  
  ### Plot re-formatted 3D structure
  r3dmol() %>%
    m_add_model(data = modified_structure_file, format = "pdb") %>%
    m_set_style(style = m_style_cartoon(
      colorfunc = paste0("
              function(atom) {
                const color = [", specific_residue_gradient,"]
                return color[Math.round(atom.b)]
              }")
    )) %>%
    m_zoom_to() 
  
}



render_structure <- function(uniprot_id,
                             file_temp_dir,
                             residue_colouring = c("Prediction quality",
                                                   "pPSE",
                                                   "Highlight all Cys",
                                                   "Specific residue(s)",
                                                   "Ligandability"),
                             residues,
                             start_col, end_col#,
                             #residue_col, rest_of_structure_col
){
  
  if(uniprot_id != ""){
    
    ## Only download if original doesn't already exist
    
    if(!dir.exists(paste0(file_temp_dir, "/", uniprot_id))){
      
      ## Download AlphaFold predicted structure
      download_alphafold_prediction(input_uniprot_id = uniprot_id,
                                    input_dir = file_temp_dir)
      
    }
    
    ## Modify downloaded structure as necessary
    # B-factor column used to plot various visualisations (prediction quality, pPSE, specific residues)
    
    if(residue_colouring == "Prediction quality"){
      
      quality_alphafold_visualise(input_uniprot_id = uniprot_id,
                                  input_dir = file_temp_dir,
                                  start_col = start_col,
                                  end_col = end_col)
      
    }else if(residue_colouring == "pPSE"){
      
      pPSE_alphafold_visualise(input_uniprot_id = uniprot_id,
                               input_dir = file_temp_dir,
                               start_col = start_col,
                               end_col = end_col)
      
    }else if(residue_colouring == "Specific residue(s)"){
      
      specific_residue_visualise(input_uniprot_id = uniprot_id,
                                 input_dir = file_temp_dir,
                                 residues,
                                 residue_col = start_col,
                                 rest_of_structure_col = end_col)
      
      
    }else if(residue_colouring == "Highlight all Cys"){
      
      all_cys_visualise(input_uniprot_id = uniprot_id,
                        input_dir = file_temp_dir,
                        residue_col = start_col,
                        rest_of_structure_col = end_col)
      
    }
    
  }
  
} # Only re-download if original structure not present




# 6. Return ligandable cysteines from given protein ID
return_ligandability_table <- function(input_uniprot_id){
  
  temp.ligand.tbl <- ligandability_df %>%
    filter(protein_id == input_uniprot_id) %>%
    as.data.frame()
  
  temp.ligand.tbl <- temp.ligand.tbl %>%
    select(protein_id, position, AA, Fragment_Name, R, quality, pPSE) %>%
    rename("UniProt ID" = 1,
           "Residue #" = 2,
           "Amino acid" = 3, 
           "Fragment" = 4,
           "Max. R Value" = 5, 
           "Prediction quality" = 6,
           "pPSE" = 7)
  
  
  return(temp.ligand.tbl)
  
  dbDisconnect()
  
}

# List of fragments we can visualise
visualise_fragments <- ligandability_df %>%
  pull(Fragment_Name) %>%
  unique() %>%
  sort() 


# 7. Print colour scale into shiny based on input colours
plot_colour_scale <- function(start_col, end_col, 
                              n_gradient){
  
  ## Use input/output colour to make a scale
  temp_gradient <- colorRampPalette(c(start_col, end_col))(n_gradient)
  
  data.frame(x = 1:n_gradient,
             y = 1,
             fill = colorRampPalette(c(start_col, 
                                       end_col))(n_gradient)) %>% 
    ggplot(aes(x = x,
               y = y,
               fill = fill,
               colour = fill)) +
    geom_tile(size=2.5) + 
    scale_fill_identity() + 
    scale_colour_identity() + 
    scale_y_continuous(expand = expansion(mult=c(0,0))) + 
    scale_x_continuous(expand = expansion(mult=c(0,0))) + 
    theme_void()
  
}

