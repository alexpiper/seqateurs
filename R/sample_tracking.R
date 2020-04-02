# Create Samplesheet ------------------------------------------------------

#' Create Samplesheet
#'
#' @param SampleSheet
#' @param runParameters
#'
#' @return
#' @export
#' @import XML
#' @import dplyr
#' @import readr
#' @import magrittr
#' @import tidyr
#' @import tibble
#' @import stringr
#' @import tibble
#'
#' @examples
create_samplesheet <- function(SampleSheets, runParameters){
  if (missing(runParameters)) {stop("Error: need to provide a runParameters file in .xml format")}
  if (missing(SampleSheets)) {stop("Error: need to provide a SampleSheet file in .csv format")}
  if (length(SampleSheets) > 1) {multi <- TRUE}
  if (!length(SampleSheets) == length(runParameters)) {stop("Error: SampleSheets and RunParameters need to be provided for every run")}

  combined = vector("list", length = length(SampleSheets))
  for (i in 1:length(SampleSheets)){

    #detect format for run
    if(any(str_detect(read_lines(runParameters[i]), "MiSeq"))){
      format <- "miseq"
      sampleskip <- 20
      header_n_max <- 19
      reads_skip <- 12
    } else if (any(str_detect(read_lines(runParameters[i]), "novaseq"))){
      format <- "novaseq"
      sampleskip = 19
      header_n_max = 18
      reads_skip = 11
    } else if (any(str_detect(read_lines(runParameters[i]), "hiseq"))){
      format <- "hiseq"
      stop("Error: HiSeq not currently supported")
    } else if (any(str_detect(read_lines(runParameters[i]), "nextseq"))){
      format <- "nextseq"
      stop("Error: NextSeq not currently supported")
    } else if (any(str_detect(read_lines(runParameters[i]), "iseq"))){
      format <- "iseq"
      stop("Error: iSeq not currently supported")
    } else(
      stop("Error: compatable platfrom not detected in runParameters file")
    )
    sample_sheet <- readr::read_csv(SampleSheets[i], skip=sampleskip, col_types = cols(
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
      sample_header <- readr::read_csv(SampleSheets[i], n_max=header_n_max, col_types = cols_only(
        `[Header]` = col_character(),
        `X2` = col_character()
      )) %>%
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

      reads <- readr::read_csv(SampleSheets[i], skip=reads_skip, n_max=2, col_types = cols_only(
        `[Reads]` = col_number() )) %>%
        pull(`[Reads]`)
      reads <- tibble::tibble(Fread = reads[1], Rread = reads[2])
    },
    warning=function(w) {if (startsWith(conditionMessage(w), "Missing column names"))
      invokeRestart("muffleWarning")})

    # Read runparameters xml
    xmlFromRunParameters <- XML::xmlParse(runParameters[i])
    run_params <- XML::xmlToDataFrame(nodes = XML::getNodeSet(xmlFromRunParameters, "/RunParameters")) %>%
      as.data.frame(stringsAsFactors=FALSE )

    if(format == "miseq"){
      run_params <- run_params %>%
        dplyr::mutate(FlowCellExpiry = FlowcellRFIDTag %>%
                        str_replace("^.{0,23}", "") %>%
                        str_replace(".{0,9}$", "") %>%
                        as.Date(),
                      ReagentKitExpiry = ReagentKitRFIDTag %>%
                        str_replace("^.{0,23}", "") %>%
                        str_replace(".{0,9}$", "") %>%
                        as.Date(),
                      PR2Expiry = PR2BottleRFIDTag %>%
                        str_replace("^.{0,23}", "") %>%
                        str_replace(".{0,9}$", "") %>%
                        as.Date(),
                      FCID = Barcode %>%
                        str_replace("^.{0,10}", ""),
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

    combined[[i]] <- sample_sheet %>%
      cbind(sample_header) %>%
      cbind(reads) %>%
      cbind(run_params)
    message("Combined sample sheets for: ")
    message(paste0(unique(combined[[i]]$FCID)," ", format, "\n"))
  }
  out <- dplyr::bind_rows(combined)
  message(paste0(length(unique(out$Sample_ID))," samples total"))
  return(out)
}
