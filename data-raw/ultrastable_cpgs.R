## code to prepare `ultrastable_cpgs`

ultrastable_cpgs = read.table("https://static-content.springer.com/esm/art%3A10.1186%2F1756-8935-7-28/MediaObjects/13072_2014_333_MOESM2_ESM.txt") |>
  tibble::rownames_to_column("probe_id") |>
  dplyr::pull(probe_id)

use_data(ultrastable_cpgs)

