# Create Samplesheet ------------------------------------------------------

#' Create Samplesheet
#'
#' @param SampleSheet A samplesheet or set of samplesheets in MiSeq or NovaSeq format
#' @param runParameters A runParamers.xml or set of these output by the sequencer
#' @param template The output format you would like.
#' This can be a character referring to the version of the sample sheet (currently only "V4" is supported)
#' Or a data.frame input to use as a template
#'
#' @return
#' @export
#' @import purrr
#' @import dplyr
#' @import janitor
#'
#' @examples
create_samplesheet <- function(SampleSheet, runParameters, template = "V4"){
  if (missing(runParameters)) {stop("Error: need to provide a runParameters file in .xml format")}
  if (missing(SampleSheet)) {stop("Error: need to provide a SampleSheet file in .csv format")}
  if (length(SampleSheet) > 1) {multi <- TRUE}
  if (!length(SampleSheet) == length(runParameters)) {
    stop("Error: you have provided ", length(SampleSheet) , " SampleSheets and ", length(runParameters), " runParameters files. One of each must be provided per run")
    }

  #Parse files
  merged <- purrr::map2(SampleSheet, runParameters, parse_seqrun) %>%
    dplyr::bind_rows()

  # Reformat to the format required
  if (is.character(template) && template=="V4"){
    # read in the template from package data
     data("samdf_template_v4", package="seqateurs")
    template <- get("samdf_template_v4", envir = sys.frame(sys.parent(0)))
  } else if (any(class(template) == "data.frame")){
    template <- template
  } else {
    stop("Error, only template='V4' or a user provided data framecurrently supported")
  }
    matching <- merged %>%
      janitor::clean_names()%>%
      dplyr::rename(
        i7_index = index,
        i5_index = index2,
        index_plate = sample_plate,
        index_well = sample_well,
        operator_name = investigator_name,
        client_name = project_name,
        seq_id = instrument_name,
        seq_date = run_start_date,
        seq_run_id = run_id
      ) %>%
      #dplyr::mutate(seq_platform = format) %>%
      dplyr::select_if(names(.) %in% colnames(template))
    matching[,setdiff(colnames(template), colnames(matching))] <- NA
    out <- matching %>%
      dplyr::select(colnames(template))

  message(paste0(length(unique(out$sample_id))," samples total"))
  return(out)
}

# Create logsheet ------------------------------------------------------

#' Create logsheet
#'
#' @param SampleSheet A samplesheet or set of samplesheets in MiSeq or NovaSeq format
#' @param runParameters A runParamers.xml or set of these output by the sequencer
#'
#' @return
#' @export
#' @import purrr
#' @import dplyr
#' @import janitor
#'
#' @examples
create_logsheet <- function(SampleSheet, runParameters){
  if (missing(runParameters)) {stop("Error: need to provide a runParameters file in .xml format")}
  if (missing(SampleSheet)) {stop("Error: need to provide a SampleSheet file in .csv format")}
  if (length(SampleSheet) > 1) {multi <- TRUE}
  if (!length(SampleSheet) == length(runParameters)) {stop("Error: SampleSheet and RunParameters need to be provided for every run")}

  #Parse files
  merged <- purrr::map2(SampleSheet, runParameters, parse_seqrun) %>%
    dplyr::bind_rows()

  out <- merged %>%
    janitor::clean_names()%>%
    dplyr::rename(
      i7_index = index,
      i5_index = index2,
      index_plate = sample_plate,
      index_well = sample_well,
      operator_name = investigator_name,
      client_name = project_name,
      seq_id = instrument_name,
      seq_date = run_start_date,
      seq_run_id = run_id
    )
  message(paste0(length(unique(out$sample_id))," samples total"))
  return(out)
}



# Parse sequencing data ---------------------------------------------------

#' Parse seqrun
#'
#' @param SampleSheet
#' @param runParameters
#'
#' @import XML
#' @import dplyr
#' @import readr
#' @import magrittr
#' @import tidyr
#' @import tibble
#' @import stringr
#' @import tibble
#' @import janitor
#' @return
#'
#' @examples
parse_seqrun <- function(SampleSheet, runParameters){
  if (missing(runParameters)) {stop("Error: need to provide a runParameters file in .xml format")}
  if (missing(SampleSheet)) {stop("Error: need to provide a SampleSheet file in .csv format")}
  if (!length(SampleSheet) == length(runParameters)) {stop("Error: SampleSheet and RunParameters need to be provided for every run")}
  #detect format for run
  if(any(stringr::str_detect(readr::read_lines(runParameters), "MiSeq"))){
    format <- "miseq"
    sampleskip <- 20
    header_n_max <- 19
    reads_skip <- 12
  } else if (any(stringr::str_detect(readr::read_lines(runParameters), "novaseq"))){
    format <- "novaseq"
    sampleskip = 19
    header_n_max = 18
    reads_skip = 11
  } else if (any(stringr::str_detect(readr::read_lines(runParameters), "hiseq"))){
    format <- "hiseq"
    stop("Error: HiSeq not currently supported")
  } else if (any(stringr::str_detect(readr::read_lines(runParameters), "nextseq"))){
    format <- "nextseq"
    stop("Error: NextSeq not currently supported")
  } else if (any(stringr::str_detect(readr::read_lines(runParameters), "iseq"))){
    format <- "iseq"
    stop("Error: iSeq not currently supported")
  } else(
    stop("Error: compatable platfrom not detected in runParameters file")
  )
  # Read in samplesheet from run
  sample_sheet <- readr::read_csv(SampleSheet, skip=sampleskip, col_types = cols(
    Sample_ID = col_character(),
    Sample_Name = col_character(),
    Sample_Plate = col_character(),
    Sample_Well = col_character(),
    I7_Index_ID = col_character(),
    index = col_character(),
    I5_Index_ID = col_character(),
    index2 = col_character(),
    Sample_Project = col_character()
  ))

  withCallingHandlers({ # Handle Annoying missing columns function
    sample_header <- readr::read_csv(SampleSheet, n_max=header_n_max) %>%
      dplyr::select(1:2) %>%
      magrittr::set_colnames(c("var", "value")) %>%
      tidyr::drop_na(var) %>%
      dplyr::mutate(var = var %>%
                      str_replace_all("InvestigatorName", "Investigator_Name") %>% #Convert camel to snake case
                      str_replace_all("ExperimentName", "Experiment_Name") %>%
                      str_replace(" ", "_") %>%
                      make.unique() %>%
                      str_replace("\\.1", "_R")
      ) %>%
      tibble::column_to_rownames("var") %>%
      t() %>%
      tibble::as_tibble() %>%
      dplyr::select_if(names(.) %in% c('Investigator_Name', 'Project_Name', 'Experiment_Name', 'Assay', 'Adapter'))

    reads <- readr::read_csv(SampleSheet, skip=reads_skip, n_max=2, col_types = cols_only(
      `[Reads]` = col_number() )) %>%
      pull(`[Reads]`)
    reads <- tibble::tibble(for_read_length = reads[1], rev_read_length = reads[2])
  },
  warning=function(w) {if (startsWith(conditionMessage(w), "Missing column names"))
    invokeRestart("muffleWarning")})

  # Read runparameters xml
  xmlFromRunParameters <- XML::xmlParse(runParameters)
  run_params <- XML::xmlToDataFrame(nodes = XML::getNodeSet(xmlFromRunParameters, "/RunParameters")) %>%
    as.data.frame(stringsAsFactors=FALSE )

  if(format == "miseq"){
    run_params <- run_params %>%
      dplyr::mutate(FlowCellExpiry = FlowcellRFIDTag %>%
                      stringr::str_replace("^.{0,23}", "") %>%
                      stringr::str_replace(".{0,9}$", "") %>%
                      as.Date(),
                    ReagentKitExpiry = ReagentKitRFIDTag %>%
                      stringr::str_replace("^.{0,23}", "") %>%
                      stringr::str_replace(".{0,9}$", "") %>%
                      as.Date(),
                    PR2Expiry = PR2BottleRFIDTag %>%
                      stringr::str_replace("^.{0,23}", "") %>%
                      stringr::str_replace(".{0,9}$", "") %>%
                      as.Date(),
                    FCID = Barcode %>%
                      stringr::str_replace("^.{0,10}", ""),
                    RunStartDate = lubridate::ymd(RunStartDate)
      ) %>%
      dplyr::rename( InstrumentName = ScannerID
      ) %>%
      dplyr::select(
        RunID,
        InstrumentName,
        RunNumber,
        FCID,
        RunStartDate,
        PR2BottleBarcode,
        ReagentKitBarcode,
        FlowCellExpiry,
        ReagentKitExpiry,
        PR2Expiry,
        MostRecentWashType) %>%
      dplyr::mutate_if(is.factor, as.character)

  } else if(format == "novaseq"){
    run_params <- run_params %>%
      dplyr::mutate(RunStartDate = lubridate::ymd(RunStartDate),
                    RunID = RunId
      ) %>%
      dplyr::select(
        RunID,
        InstrumentName,
        RunNumber,
        RunStartDate) %>%
      dplyr::mutate_if(is.factor, as.character)

    RFIDS <- XML::xmlToDataFrame(nodes = XML::getNodeSet(xmlFromRunParameters, "//RfidsInfo"))%>%
      as.data.frame(stringsAsFactors=FALSE) %>%
      dplyr::rename(
        FCID = FlowCellSerialBarcode,
        LibTubeID = LibraryTubeSerialBarcode,
        SbsID = SbsSerialBarcode,
        ClusterID = ClusterSerialBarcode,
        BufferID = BufferSerialBarcode,
      ) %>%
      dplyr::mutate(
        FlowCellExpiry = lubridate::mdy(str_remove(FlowCellExpirationdate," 00:00:00")),
        SbsExpiry = lubridate::mdy(str_remove(SbsExpirationdate," 00:00:00")),
        ClusterExpiry = lubridate::mdy(str_remove(ClusterExpirationdate," 00:00:00")),
        BufferExpiry = lubridate::mdy(str_remove(BufferExpirationdate," 00:00:00")),
      )%>%
      dplyr::select(
        FCID,
        LibTubeID,
        ClusterID,
        SbsID,
        BufferID,
        FlowCellExpiry,
        SbsExpiry,
        ClusterExpiry,
        BufferExpiry
      ) %>%
      dplyr::mutate_if(is.factor, as.character)
    run_params <- dplyr::bind_cols(run_params, RFIDS)
  }

  #Merge different sheets
  combined <- sample_sheet %>%
    cbind(sample_header) %>%
    cbind(reads) %>%
    cbind(run_params)
  message("Combined sample sheets for: ")
  message(paste0(unique(combined[]$FCID)," ", format, "\n"))
  return(combined)
}
