## Read in packages, data (in form of sqlite DBs and one dataframe) and functions for plotting/structure manipulation
source("load_data_and_functions.R")


# Create the UI
ui <- fluidPage(
  
  title = "AlpACA-DB",
  
  # App title ----
  fluidRow(style = "background-color:#ECECEC;",
           column(12,
                  h1(strong("Alp", .noWS = c("before", "after")), "haFold2 ", 
                     strong("A", .noWS = c("before", "after")), "ccessibility and ", 
                     strong("C", .noWS = c("before", "after")), "ys ", 
                     strong("A", .noWS = c("before", "after")), "nnotation ", 
                     strong("D", .noWS = c("before", "after")), "ata", 
                     strong("b", .noWS = c("before", "after")), "ase"), .noWS = c("before", "after"))),
  
  fluidRow(
    style = "background-color:#ECECEC;",
    column(9, 
           h4("Tool for exploring ligandable cysteines on >20,000 human predicted protein structures, accompanying an", 
              a("accessibility-based analysis", href="https://www.biorxiv.org/content/10.1101/2022.12.12.518491v1"), 
              " of published cysteine-focused chemoproteomic datasets. "),
           br())),
  
  tags$head(tags$style(
    type="text/css",
    "#fragment_structure img {max-width: 75%; width: 75%; max-height: 50%; height: auto}"
  )),
  
  
  fluidRow(
    style="padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:10px;background-color:#ECECEC;",
    column(3,
           h4("UniProt ID"),
           textInput("uniprotID", 
                     NULL),
           span(textOutput("uniprot_input_infotext_1"), 
                style="color: grey; font-size:13px; opacity: 0.75"),
           br(),
           span(textOutput("uniprot_input_infotext_2"), 
                style="color: grey; font-size:13px; opacity: 0.75"),
           br(),
           br()),
    
    column(3,
           h4("Visualisation type"),
           selectInput("select_vis", 
                       NULL, 
                       choices = c("Select visualisation",
                                   "pPSE",
                                   "Prediction quality",
                                   "Highlight all Cys",
                                   "Specific residue(s)"),
                       selected = "Select visualisation"),
           #uiOutput("conditional_pPSE_slider"),
           uiOutput("specific_residue"),
           span(textOutput("specific_residue_info"), 
                style="color: grey; font-size:13px; opacity: 0.75"),
           span(textOutput("pPSE_info"), 
                style="color: grey; font-size:13px; opacity: 0.75"),
           span(textOutput("prediction_quality_info"), 
                style="color: grey; font-size:13px; opacity: 0.75"),
           span(textOutput("visualisation_info"), 
                style="color: grey; font-size:13px; opacity: 0.75"),
           uiOutput("conditional_ligand"),
           uiOutput("conditional_quality")),
    
    column(3,
           h4("Colour scheme"),
           uiOutput("conditional_col1"),
           uiOutput("conditional_col2"),
           plotOutput("colour_scale", width = "200px", height = "20px")),
    
    column(3,
           h4("Extra info"),
           span(textOutput("colour_info_1"), 
                style="color: grey; font-size:13px; opacity: 0.75"),
           h3("    "),
           span(textOutput("colour_info_2"), 
                style="color: grey; font-size:13px; opacity: 0.75")
           
    )
    
    
  ),
  
  fluidRow(style="padding-left:5px; padding-right:5px; padding-top:30px; padding-bottom:15px",
           column(9,
                  r3dmolOutput("structurePlot")),
           column(3, 
                  style = "background-color: white;",
                  h4("Fragment structure"),
                  selectInput("frag_name",
                              NULL,
                              choices = c("Select fragment",
                                          visualise_fragments),
                              selected = "Select fragment"),
                  imageOutput("fragment_structure"))),
  
  fluidRow(style="padding-left:5px; padding-right:5px; padding-top:10px; padding-bottom:10px",
           column(12,
                  shiny::dataTableOutput("ligandability_table"))),
  
  
  fluidRow(style="padding-left:15px",
           h5("Tate Lab (2023). Please leave any issues/feedback on our ",
              a("GitHub", href="https://github.com/TateLab/"),
              "repository."))
  
)

# Create the server
server <- function(input, output) {
  
  output$uniprot_input_infotext_1 <- renderText({
    
    "1. Enter a UniProt accession number."
    
  })
  
  output$uniprot_input_infotext_2 <- renderText({
    
    "Note, currently only human proteins in the AlphaFold v4 database are supported."
    
  })
  
  output$ligandability_table <- shiny::renderDataTable({
    
    validate(need(nrow(return_ligandability_table(input$uniprotID)) != 0, "No ligandability data available"))
    
    return_ligandability_table(input$uniprotID)
    
  },
  
  options = list(scrollX = TRUE,
                 scrollY = "250px")
  
  )
  
  ### Render image of fragment structure
  output$fragment_structure <- renderImage({
    req(input$frag_name != "Select fragment")
    
    file_location <- dir_ls("shiny_data/fragment_structures/", regexp = paste0(input$frag_name, ".png"))
    
    list(src = file_location,  
         contentType = 'image/png')},
    
    deleteFile = F)
  
  ### Render 3D structure
  output$structurePlot <- renderR3dmol({
    validate(
      need(RCurl::url.exists(paste0("https://alphafold.ebi.ac.uk/files/AF-", 
                                    input$uniprotID,
                                    "-F1-model_v4.pdb")), "Enter a valid UniProt ID."),
      need(input$select_vis != "Select visualisation", "Select visualisation type from drop down menu")
    )
    
    ## Conditional on input
    expr = render_structure(uniprot_id = input$uniprotID,
                            file_temp_dir = tempdir(),
                            residue_colouring = input$select_vis,
                            start_col = input$choose_col1,
                            end_col = input$choose_col2,
                            residues = input$input_specific_residue)
  })
  
  
  
  ### Output text below showing protein name
  #output$text_output <- renderPrint({
  #    
  #  })
  
  ### Conditional appearance of column 3 options
  
  output$conditional_col1 <- renderUI({
    #req(input$select_vis == "pPSE")
    
    textInput("choose_col1", 
              NULL,
              value = "blue")
    
    
  })
  
  output$conditional_col2 <- renderUI({
    #req(input$select_vis == "pPSE")
    
    textInput("choose_col2", 
              NULL, 
              value = "white")
    
  })
  
  
  
  # Ligandability
  
  output$colour_scale <- renderPlot({
    validate(
      need(!is.na(input$choose_col1), "Test1"),
      need(!is.na(input$choose_col1), "Test2")
    )
    
    plot_colour_scale(start_col = input$choose_col1,
                      end_col = input$choose_col2,
                      n_gradient = 101)
    
  })
  
  # Quality
  #output$conditional_quality <- renderUI({
  #  req(input$select_vis == "Prediction quality")
  #  
  #  sliderInput("slider1", h5("Prediction quality"),
  #              min = 20, max = 100, value = 1)
  #})
  
  # Specific residue
  output$specific_residue <- renderUI({
    req(input$select_vis == "Specific residue(s)")
    
    textInput("input_specific_residue", 
              label = NULL,
              value = "Residue #")
    
  })
  
  output$specific_residue_info <- renderText({
    req(input$select_vis == "Specific residue(s)")
    
    "Provide residue number(s) from UniProt sequence, either individual numbers, separated by commas (e.g. '152, 156') or as a range (e.g. '1-50')"
    
  })
  
  output$pPSE_info <- renderText({
    req(input$select_vis == "pPSE")
    
    "Prediction-aware part-sphere exposure (pPSE) as calculated in the corresponding pre-print"
    
  })
  
  output$prediction_quality_info <- renderText({
    req(input$select_vis == "Prediction quality")
    
    "Prediction quality (pLDDT) per residue from AlphaFold structure prediction"
    
  })
  
  output$visualisation_info <- renderText({
    req(input$select_vis == "Select visualisation")
    
    "2. Choose a visualisation to display AlphaFold structure"
    
  })
  
  output$colour_info_1 <- renderText({
    
    "2-colour scale can be customised by either specifying specific colour or hex codes (e.g. 'red', '#FFFFFF')."
    
  })
  
  output$colour_info_2 <- renderText({
    
    "Drag and scroll with mouse to manipulate the displayed structure."
    
  })
  
}

# Create the Shiny app
shinyApp(ui = ui, server = server)

