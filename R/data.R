#' Bodyfat
#'
#' The response (`y`) corresponds to
#' estimates of percentage of body fat from application of
#' Siri's 1956 equation to measurements of underwater weighing, as well as
#' age, weight, height, and a variety of
#' body circumference measurements.
#'
#' @format A list with two items representing 252 observations from
#'   14 variables
#' \describe{
#'   \item{age}{age (years)}
#'   \item{weight}{weight (lbs)}
#'   \item{height}{height (inches)}
#'   \item{neck}{neck circumference (cm)}
#'   \item{chest}{chest circumference (cm)}
#'   \item{abdomen}{abdomen circumference (cm)}
#'   \item{hip}{hip circumference (cm)}
#'   \item{thigh}{thigh circumference (cm)}
#'   \item{knee}{knee circumference (cm)}
#'   \item{ankle}{ankle circumference (cm)}
#'   \item{biceps}{biceps circumference (cm)}
#'   \item{forearm}{forearm circumference (cm)}
#'   \item{wrist}{wrist circumference (cm)}
#' }
#' @source http://lib.stat.cmu.edu/datasets/bodyfat
#' @source https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression.html
#' @family datasets
"bodyfat"

#' Abalone
#'
#' This data set contains observations of abalones, the common
#' name for any of a group of sea snails. The goal is to predict the
#' age of an individual abalone given physical measurements such as
#' sex, weight, and height.
#'
#' Only a stratified sample of 211 rows of the original data set are used here.
#'
#' @format A list with two items representing 211 observations from
#'   9 variables
#' \describe{
#'   \item{sex}{sex of abalone, 1 for female}
#'   \item{infant}{indicates that the person is an infant}
#'   \item{length}{longest shell measurement in mm}
#'   \item{diameter}{perpendicular to length in mm}
#'   \item{height}{height in mm including meat in shell}
#'   \item{weight_whole}{weight of entire abalone}
#'   \item{weight_shucked}{weight of meat}
#'   \item{weight_viscera}{weight of viscera}
#'   \item{weight_shell}{weight of shell}
#'   \item{rings}{rings. +1.5 gives the age in years}
#' }
#' @source Pace, R. Kelley and Ronald Barry, Sparse Spatial Autoregressions,
#'   Statistics and Probability Letters, 33 (1997) 291-297.
#' @family datasets
"abalone"

#' Heart disease
#'
#' Diagnostic attributes of patients classified as having heart disease or not.
#'
#' @section Preprocessing:
#' The original dataset contained 13 variables. The nominal of these were
#' dummycoded, removing the first category. No precise information regarding
#' variables `chest_pain`, `thal` and `ecg` could be found, which explains
#' their obscure definitions here.
#'
#' @format 270 observations from 17 variables represented as a list consisting
#' of a binary factor response vector `y`,
#' with levels 'absence' and 'presence' indicating the absence or presence of
#' heart disease and `x`: a sparse feature matrix of class 'dgCMatrix' with the
#' following variables:
#' \describe{
#'   \item{age}{age}
#'   \item{bp}{diastolic blood pressure}
#'   \item{chol}{serum cholesterol in mg/dl}
#'   \item{hr}{maximum heart rate achieved}
#'   \item{old_peak}{ST depression induced by exercise relative to rest}
#'   \item{vessels}{the number of major blood vessels (0 to 3) that were
#'                  colored by fluoroscopy}
#'   \item{sex}{sex of the participant: 0 for male, 1 for female}
#'   \item{angina}{a dummy variable indicating whether the person suffered
#'                 angina-pectoris during exercise}
#'   \item{glucose_high}{indicates a fasting blood sugar over 120 mg/dl}
#'   \item{cp_typical}{typical angina}
#'   \item{cp_atypical}{atypical angina}
#'   \item{cp_nonanginal}{non-anginal pain}
#'   \item{ecg_abnormal}{indicates a ST-T wave abnormality
#'                       (T wave inversions and/or ST elevation or depression of
#'                       > 0.05 mV)}
#'   \item{ecg_estes}{probable or definite left ventricular hypertrophy by
#'                    Estes' criteria}
#'   \item{slope_flat}{a flat ST curve during peak exercise}
#'   \item{slope_downsloping}{a downwards-sloping ST curve during peak exercise}
#'   \item{thal_reversible}{reversible defect}
#'   \item{thal_fixed}{fixed defect}
#' }
#' @source Dua, D. and Karra Taniskidou, E. (2017). UCI Machine Learning Repository
#'   <http://archive.ics.uci.edu/ml/>. Irvine, CA: University of California,
#'   School of Information and Computer Science.
#' @source <https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary.html#heart>
#' @family datasets
"heart"

#' Wine cultivars
#'
#' A data set of results from chemical analysis of wines grown in Italy
#' from three different cultivars.
#'
#' @format 178 observations from 13 variables represented as a list consisting
#' of a categorical response vector `y`
#' with three levels: *A*, *B*, and *C* representing different
#' cultivars of wine as well as `x`: a sparse feature matrix of class
#' 'dgCMatrix' with the following variables:
#' \describe{
#'   \item{alcohol}{alcoholic content}
#'   \item{malic}{malic acid}
#'   \item{ash}{ash}
#'   \item{alcalinity}{alcalinity of ash}
#'   \item{magnesium}{magnemium}
#'   \item{phenols}{total phenols}
#'   \item{flavanoids}{flavanoids}
#'   \item{nonflavanoids}{nonflavanoid phenols}
#'   \item{proanthocyanins}{proanthocyanins}
#'   \item{color}{color intensity}
#'   \item{hue}{hue}
#'   \item{dilution}{OD280/OD315 of diluted wines}
#'   \item{proline}{proline}
#' }
#'
#' @source Dua, D. and Karra Taniskidou, E. (2017). UCI Machine Learning Repository
#'   <http://archive.ics.uci.edu/ml/>. Irvine, CA: University of California,
#'   School of Information and Computer Science.
#' @source <https://raw.githubusercontent.com/hadley/rminds/master/1-data/wine.csv>
#' @source <https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multiclass.html#wine>
#' @family datasets
"wine"

#' Student performance
#'
#' A data set of the attributes of 382 students in secondary education
#' collected from two schools. The goal is to predict the
#' grade in math and Portugese at the end of the third period. See the
#' cited sources for additional information.
#'
#' @section Preprocessing:
#' All of the grade-specific predictors were dropped from the data set.
#' (Note that it is not clear from the source why some of these predictors are
#' specific to each grade, such as which parent is the student's guardian.)
#' The categorical variables were dummy-coded. Only the final grades (G3)
#' were kept as dependent variables, whilst the
#' first and second period grades were dropped.
#'
#' @format 382 observations from 13 variables represented as a list consisting
#' of a binary factor response matrix `y` with two responses: `portugese` and
#' `math` for the final scores in period three for the respective subjects.
#' The list also contains `x`: a sparse feature matrix of class
#' 'dgCMatrix' with the following variables:
#' \describe{
#'   \item{school_ms}{student's primary school, 1 for Mousinho da Silveira and 0
#'         for Gabriel Pereira}
#'   \item{sex}{sex of student, 1 for male}
#'   \item{age}{age of student}
#'   \item{urban}{urban (1) or rural (0) home address}
#'   \item{large_family}{whether the family size is larger than 3}
#'   \item{cohabitation}{whether parents live together}
#'   \item{Medu}{mother's level of education (ordered)}
#'   \item{Fedu}{fathers's level of education (ordered)}
#'   \item{Mjob_health}{whether the mother was employed in health care}
#'   \item{Mjob_other}{whether the mother was employed as something other than
#'         the specified job roles}
#'   \item{Mjob_services}{whether the mother was employed in the service sector}
#'   \item{Mjob_teacher}{whether the mother was employed as a teacher}
#'   \item{Fjob_health}{whether the father was employed in health care}
#'   \item{Fjob_other}{whether the father was employed as something other than
#'         the specified job roles}
#'   \item{Fjob_services}{whether the father was employed in the service sector}
#'   \item{Fjob_teacher}{whether the father was employed as a teacher}
#'   \item{reason_home}{school chosen for being close to home}
#'   \item{reason_other}{school chosen for another reason}
#'   \item{reason_rep}{school chosen for its reputation}
#'   \item{nursery}{whether the student attended nursery school}
#'   \item{internet}{Pwhether the student has internet access at home}
#' }
#'
#' @source P. Cortez and A. Silva. Using Data Mining to Predict Secondary School
#'   Student Performance. In A. Brito and J. Teixeira Eds., Proceedings of 5th
#'   FUture BUsiness TEChnology Conference (FUBUTEC 2008) pp. 5-12, Porto,
#'   Portugal, April, 2008, EUROSIS, ISBN 978-9077381-39-7.
#' @source Dua, D. and Karra Taniskidou, E. (2017). UCI Machine Learning
#'  Repository <http://archive.ics.uci.edu/ml/>. Irvine, CA: University of
#'  California, School of Information and Computer Science.
#' @family datasets
"student"

#' Glioma metabolomics
#'
#' Metabolomics dataset from 165 different plasma measurements from
#' 94 patients (cases) with glioma (brain tumours) and 71 healthy controls. The
#' goal is to predict whether a sample is from a patient or a control based
#' on the metabolite measurements.
#'
#' @section Preprocessing:
#' We have removed the patients with meningioma from the original dataset (which
#' contained 235 samples) to create a binary classification problem. Also,
#' the authors originally had 188 features but removed some of these due to
#' missing data.
#'
#' @format 165 observations from 138 variables represented as a list consisting
#'   of a binary response (factor) vector `y` with levels 'control' and 'case'
#'   indicating whether the sample is from a healthy control or a patient with
#'   glioma, as well as `x`: a matrix of 138 metabolite measurements.
#'
#' @source Godlewski, A., Czajkowski, M., Mojsak, P., Pienkowski, T., Gosk, W.,
#'   Lyson, T., Mariak, Z., Reszec, J., Kondraciuk, M., Kaminski, K., Kretowski,
#'   M., Moniuszko, M., Kretowski, A., & Ciborowski, M. (2023). A comparison of
#'   different machine-learning techniques for the selection of a panel of
#'   metabolites allowing early detection of brain tumors. Scientific Reports,
#'   13(1), 11044. \doi{10.1038/s41598-023-38243-1}
"glioma"
