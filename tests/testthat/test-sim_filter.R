test_that("`sim_filter_keep_or_drop_some` works", {
  n <- 5

  set.seed(42)
  population <- tibble::tibble(
    Metadata_group = sample(c("a", "b"), n, replace = TRUE),
    Metadata_type = sample(c("x", "y"), n, replace = TRUE),
    x = rnorm(n),
    y = x + rnorm(n) / 100,
    z = y + rnorm(n) / 1000
  )
  annotation_cols <- c("Metadata_group", "Metadata_type")
  # this is not a great test because it tests more than one function
  sim_df <- matric::sim_calculate(population, method = "pearson")
  row_metadata <- attr(sim_df, "row_metadata")
  sim_df <- matric::sim_annotate(sim_df, row_metadata, annotation_cols)
  filter_keep <-
    tibble::tibble(Metadata_group = "a", Metadata_type = "x")
  filter_drop <-
    tibble::tibble(Metadata_group = "a", Metadata_type = "x")

  s1 <-
    matric::sim_filter_keep_or_drop_some(sim_df, row_metadata, filter_keep = filter_keep, filter_side = "left")
  s2 <-
    matric::sim_filter_keep_or_drop_some(sim_df, row_metadata, filter_drop = filter_drop, filter_side = "left")

  expect_equal(
    dplyr::bind_rows(s1, s2) %>% dplyr::arrange(id1, id2),
    sim_df %>% dplyr::arrange(id1, id2)
  )

  s1 <-
    matric::sim_filter_keep_or_drop_some(sim_df, row_metadata, filter_keep = NULL, filter_side = "left")
  s2 <-
    matric::sim_filter_keep_or_drop_some(sim_df, row_metadata, filter_drop = NULL, filter_side = "left")


  expect_equal(
    s1, s2
  )
})

test_that("`sim_filter_all_same` works", {
  n <- 5

  set.seed(42)
  population <- tibble::tibble(
    Metadata_group = sample(c("a", "b"), n, replace = TRUE),
    Metadata_type = sample(c("x", "y"), n, replace = TRUE),
    x = rnorm(n),
    y = x + rnorm(n) / 100,
    z = y + rnorm(n) / 1000
  )
  annotation_cols <- c("Metadata_group", "Metadata_type")

  sim_df <- matric::sim_calculate(population, method = "pearson")
  row_metadata <- attr(sim_df, "row_metadata")
  sim_df <- matric::sim_annotate(sim_df, row_metadata, annotation_cols)

  all_same_cols <- c("Metadata_group")
  include_group_tag <- TRUE
  drop_lower <- FALSE

  sim_df <- matric::sim_filter_all_same(
    sim_df,
    row_metadata,
    all_same_cols,
    annotation_cols,
    include_group_tag,
    drop_lower
  )

  sim_df <-
    matric::sim_annotate(sim_df, row_metadata, annotation_cols = annotation_cols)

  expect_equal(sim_df$Metadata_group1, sim_df$Metadata_group2)
})


test_that("`sim_filter_all_same_keep_some` works", {
  n <- 20

  set.seed(42)
  population <- tibble::tibble(
    Metadata_group = sample(c("a", "b"), n, replace = TRUE),
    Metadata_type = sample(c("x", "y"), n, replace = TRUE),
    x = rnorm(n),
    y = x + rnorm(n) / 100,
    z = y + rnorm(n) / 1000
  )
  annotation_cols <- c("Metadata_group", "Metadata_type")

  sim_df <- matric::sim_calculate(population, method = "pearson")
  row_metadata <- attr(sim_df, "row_metadata")
  sim_df <- matric::sim_annotate(sim_df, row_metadata, annotation_cols)

  all_same_cols <- c("Metadata_group")
  filter_keep_right <-
    tibble::tibble(Metadata_group = "a", Metadata_type = "x")
  drop_reference <- FALSE

  sim_df <- matric::sim_filter_all_same_keep_some(
    sim_df,
    row_metadata,
    all_same_cols,
    filter_keep_right,
    annotation_cols,
    drop_reference
  )

  sim_df <-
    matric::sim_annotate(sim_df, row_metadata, annotation_cols = annotation_cols)

  expect_equal(sim_df$Metadata_group1, sim_df$Metadata_group2)

  expect_equal(
    sim_df %>%
      dplyr::distinct(Metadata_group2, Metadata_type2) %>%
      as.character(),
    filter_keep_right %>%
      as.character()
  )
})


test_that("sim_filter_some_different_drop_some works", {
  n <- 20

  set.seed(42)
  population <- tibble::tibble(
    Metadata_group = sample(c("a", "b"), n, replace = TRUE),
    Metadata_type1 = sample(c("x", "y"), n, replace = TRUE),
    Metadata_type2 = sample(c("p", "q"), n, replace = TRUE),
    Metadata_type3 = sample(c("r", "s"), n, replace = TRUE),
    x = rnorm(n),
    y = x + rnorm(n) / 100,
    z = y + rnorm(n) / 1000
  )
  annotation_cols <-
    c(
      "Metadata_group",
      "Metadata_type1",
      "Metadata_type2",
      "Metadata_type3"
    )

  sim_df <- matric::sim_calculate(population, method = "pearson")
  row_metadata <- attr(sim_df, "row_metadata")
  sim_df <- matric::sim_annotate(sim_df, row_metadata, annotation_cols)

  all_same_cols <- c("Metadata_group")
  all_different_cols <- c("Metadata_type1")
  any_different_cols <- c("Metadata_type2", "Metadata_type3")

  filter_drop_left <-
    tibble::tibble(Metadata_group = "a", Metadata_type1 = "x")
  filter_drop_right <-
    tibble::tibble(Metadata_group = "a", Metadata_type1 = "x")
  drop_reference <- FALSE

  sim_df <-
    matric::sim_filter_some_different_drop_some(
      sim_df,
      row_metadata,
      any_different_cols,
      all_same_cols,
      all_different_cols,
      filter_drop_left,
      filter_drop_right,
      annotation_cols
    )

  sim_df <-
    matric::sim_annotate(sim_df, row_metadata, annotation_cols = annotation_cols)

  expect_equal(sim_df$Metadata_group1, sim_df$Metadata_group2)

  expect_true(all(sim_df$Metadata_type11 != sim_df$Metadata_type12))

  expect_true(all(
    (sim_df$Metadata_type21 != sim_df$Metadata_type22) |
      (sim_df$Metadata_type31 != sim_df$Metadata_type32)
  ))
})
