#' Two channels of ECG measurments obtained from the QT database
#'
#'
#' Data from record sel102. 45,000 measurments of 2 ECG channels over 180 seconds.
#'
#' @name ECG_measurements
#'
#' @docType data
#'
#' @usage data(ECG_measurements)
#'
#' @format A data frame with 45000 rows and 3 variables:
#' \describe{
#'   \item{time}{Time of measurment in seconds}
#'   \item{channel1}{ECG reading of first channel}
#'   \item{channel2}{ECG reading of second channel}
#'   }
#'
#' @keywords datasets
#'
#' @references Laguna, P., Mark, R. G., Goldberg, A., & Moody, G. B.
#'  (1997, September). A database for evaluation of algorithms for measurement
#'   of QT and other waveform intervals in the ECG. In Computers in cardiology
#'    1997 (pp. 673-676). IEEE.
#'
#' @source \href{https://physionet.org/content/qtdb/1.0.0/}{QT database}
#'
#' @examples
#' \donttest{
#' data(ECG_measurements)
#' }
NULL
