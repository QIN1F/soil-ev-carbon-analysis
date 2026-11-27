# ============================================================
# MAG 所代表分类单元的全球塑料际分布图（完整世界地图 v3）
# 分析层级：
#   1) Strict Genus：仅允许本地 Genus 与全球 Genus 完全同名
#   2) Family：仅允许本地 Family 与全球 Family 完全同名
#
# 地图风格参考 Qin et al. (2026) One Earth, Figure 3E：
# 灰色世界底图 + 半透明暖色气泡 + 左侧图例。
#
# 数据处理原则：
#   - 删除本地 UBA5507、UBA6024 和 Cyanobacteria MAG
#   - 全球丰度表按 MAGs 去重，避免一个 MAG 因多个 Subject id 重复累加
#   - 样品相对丰度 =
#       目标分类单元对应唯一 MAG 丰度之和 /
#       样品内全部唯一 MAG 丰度之和 × 100
#   - 地图范围固定为 180°W–180°E、90°S–90°N，完整显示世界地图
#   - 地图按“唯一经纬度”聚合，不对坐标做 jitter 或人工移动
#   - 同一经纬度的重复样品数控制气泡大小
#   - 同一经纬度所有样品的平均相对丰度控制气泡颜色
#
# 输入文件应放在当前工作目录：
#   MAG-result-0(1).csv
#   bin_host_abundance_20260508.txt
#   Sample_site_group_new.txt
# ============================================================

# ---------- 0. 安装并加载软件包 ----------
# %>% 由 magrittr 提供；dplyr 会重新导出该操作符。
# 这里显式安装并加载各依赖包，避免只运行部分代码时找不到 %>%。

options(repos = c(CRAN = "https://cloud.r-project.org"))

required_packages <- c(
  "magrittr", "dplyr", "tidyr", "readr", "tibble",
  "purrr", "stringr", "ggplot2", "maps", "scales"
)

new_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(new_packages) > 0) {
  install.packages(new_packages, dependencies = TRUE)
}

# 检查安装是否成功
failed_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(failed_packages) > 0) {
  stop(
    "以下R包安装失败：",
    paste(failed_packages, collapse = ", "),
    "\n请检查网络、CRAN镜像或R版本后重新运行。"
  )
}

suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(maps)
  library(scales)
})

# 显式注册管道操作符，保证后续 local_mag %>% mutate(...) 可以运行
`%>%` <- magrittr::`%>%`

message("R包已成功加载，%>% 管道操作符可用。")

# ---------- 1. 文件路径和用户设置 ----------
mag_file <- "MAG-result.csv"
abundance_file <- "bin_host_abundance_20260508.txt"
site_file <- "Sample_site_group_new.txt"

output_dir <- "MAG_global_map_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 默认设置：气泡大小 = 同一经纬度的样品数。
# 如需完全按文献 Figure 3E 图注，让气泡大小也表示丰度，
# 改为 bubble_size_by <- "abundance"
bubble_size_by <- "sample_count"  # 可选："sample_count" 或 "abundance"

# 本地 MAG 的删除规则
excluded_genera <- c("UBA5507", "UBA6024")
excluded_phyla <- c("Cyanobacteria")

# 暖色渐变，接近参考图
bubble_palette <- c(
  "#0E6DB3", "#5CAAD7", "#95C6DE",
  "#F6BCA9", "#D35B4D", "#BB1E38", "#8E0D29"
)

# ---------- 2. 辅助函数 ----------
clean_rank <- function(x, prefix) {
  x <- trimws(as.character(x))
  x <- sub(paste0("^", prefix), "", x)
  x[x == ""] <- NA_character_
  x
}

positive_min <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) NA_real_ else min(x)
}

positive_max <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) NA_real_ else max(x)
}

collapse_unique <- function(x) {
  x <- sort(unique(x[!is.na(x) & x != ""]))
  paste(x, collapse = ", ")
}

format_lon <- function(x) {
  ifelse(
    x < 0, paste0(abs(x), "°W"),
    ifelse(x > 0, paste0(x, "°E"), "0°")
  )
}

format_lat <- function(x) {
  ifelse(
    x < 0, paste0(abs(x), "°S"),
    ifelse(x > 0, paste0(x, "°N"), "0°")
  )
}

# ---------- 3. 读取本地 MAG ----------
local_mag <- readr::read_csv(
  mag_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)

required_local_cols <- c("MAG", "Phylum", "Family", "Genus")
missing_local_cols <- setdiff(required_local_cols, names(local_mag))

if (length(missing_local_cols) > 0) {
  stop("本地 MAG 表缺少列：", paste(missing_local_cols, collapse = ", "))
}

local_mag <- local_mag %>%
  mutate(
    Phylum = clean_rank(Phylum, "p__"),
    Family = clean_rank(Family, "f__"),
    Genus = clean_rank(Genus, "g__")
  )

local_mag_keep <- local_mag %>%
  filter(
    !(Genus %in% excluded_genera),
    is.na(Phylum) | !(Phylum %in% excluded_phyla)
  )

# ---------- 4. 读取样品经纬度，并进行严格坐标核验 ----------
site_meta <- readr::read_tsv(
  site_file,
  show_col_types = FALSE,
  name_repair = "minimal"
) %>%
  mutate(
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
  )

required_site_cols <- c(
  "Sample_name", "Site", "longitude", "latitude",
  "sample_type", "Region", "Country"
)
missing_site_cols <- setdiff(required_site_cols, names(site_meta))

if (length(missing_site_cols) > 0) {
  stop("位点表缺少列：", paste(missing_site_cols, collapse = ", "))
}

if (anyDuplicated(site_meta$Sample_name) > 0) {
  duplicated_samples <- site_meta %>%
    count(Sample_name) %>%
    filter(n > 1) %>%
    pull(Sample_name)
  stop(
    "位点表中的 Sample_name 不唯一：",
    paste(duplicated_samples, collapse = ", ")
  )
}

if (any(!is.finite(site_meta$longitude)) ||
    any(!is.finite(site_meta$latitude))) {
  stop("位点表存在无法解析的经纬度。")
}

if (any(site_meta$longitude < -180 | site_meta$longitude > 180)) {
  stop("存在超出 [-180, 180] 的经度。")
}

if (any(site_meta$latitude < -90 | site_meta$latitude > 90)) {
  stop("存在超出 [-90, 90] 的纬度。")
}

# 同一个 Site ID 必须只有一个经纬度
site_coord_audit <- site_meta %>%
  group_by(Site) %>%
  summarise(
    n_longitude = n_distinct(longitude),
    n_latitude = n_distinct(latitude),
    longitude = first(longitude),
    latitude = first(latitude),
    .groups = "drop"
  )

if (any(site_coord_audit$n_longitude != 1 |
        site_coord_audit$n_latitude != 1)) {
  stop("至少一个 Site ID 对应多个不同经纬度，请先核对位点表。")
}

# 不同 Site ID 可以具有完全相同的坐标。
# 地图必须按唯一经纬度聚合，否则多个气泡会完全重叠。
location_lookup <- site_meta %>%
  distinct(longitude, latitude) %>%
  arrange(longitude, latitude) %>%
  mutate(Location_ID = sprintf("Location_%03d", row_number()))

site_meta <- site_meta %>%
  left_join(
    location_lookup,
    by = c("longitude", "latitude"),
    relationship = "many-to-one"
  )

# 同一经纬度下的样品类型、国家应一致
location_coord_audit <- site_meta %>%
  group_by(Location_ID, longitude, latitude) %>%
  summarise(
    n_sample_type = n_distinct(sample_type),
    n_region = n_distinct(Region),
    n_country = n_distinct(Country),
    Site_IDs = collapse_unique(Site),
    sample_size = n_distinct(Sample_name),
    sample_type = first(sample_type),
    Region = first(Region),
    Country = first(Country),
    .groups = "drop"
  )

if (any(location_coord_audit$n_sample_type != 1 |
        location_coord_audit$n_region != 1 |
        location_coord_audit$n_country != 1)) {
  stop("同一经纬度下存在冲突的样品类型、Region 或 Country。")
}

# 输出坐标核验表
readr::write_csv(
  location_coord_audit,
  file.path(output_dir, "coordinate_audit_unique_locations.csv")
)

# ---------- 5. 读取全球丰度表 ----------
global_raw <- readr::read_tsv(
  abundance_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)

taxonomy_cols <- c(
  "MAGs", "Subject id", "Domain", "Phylum", "Class",
  "Order", "Family", "Genus", "Species"
)

missing_global_cols <- setdiff(taxonomy_cols, names(global_raw))
if (length(missing_global_cols) > 0) {
  stop("全球丰度表缺少列：", paste(missing_global_cols, collapse = ", "))
}

sample_cols <- setdiff(names(global_raw), taxonomy_cols)

# 丰度表中的样品必须与位点表完全对应
missing_meta_samples <- setdiff(sample_cols, site_meta$Sample_name)
extra_meta_samples <- setdiff(site_meta$Sample_name, sample_cols)

if (length(missing_meta_samples) > 0) {
  stop(
    "下列丰度样品缺少经纬度：",
    paste(missing_meta_samples, collapse = ", ")
  )
}

if (length(extra_meta_samples) > 0) {
  warning(
    "位点表中存在未进入丰度表的样品：",
    paste(extra_meta_samples, collapse = ", ")
  )
}

global_raw <- global_raw %>%
  mutate(
    Phylum_clean = clean_rank(Phylum, "p__"),
    Family_clean = clean_rank(Family, "f__"),
    Genus_clean = clean_rank(Genus, "g__")
  ) %>%
  mutate(
    across(
      all_of(sample_cols),
      ~ replace_na(suppressWarnings(as.numeric(.x)), 0)
    )
  )

# ---------- 6. 检查并删除重复 MAG ----------
# 同一 MAGs 如果因多个 Subject id 重复出现，
# 分类和全部样品丰度必须完全相同，才允许保留第一条。
profile_matrix <- global_raw %>%
  select(all_of(sample_cols)) %>%
  as.matrix()

profile_key <- apply(
  profile_matrix,
  1,
  function(x) paste(format(x, digits = 15, scientific = FALSE, trim = TRUE),
                     collapse = "|")
)

duplicate_check <- global_raw %>%
  mutate(
    .taxonomy_key = paste(
      Phylum_clean, Family_clean, Genus_clean, Species,
      sep = "|"
    ),
    .profile_key = profile_key
  ) %>%
  group_by(MAGs) %>%
  summarise(
    n_rows = n(),
    n_taxonomy_profiles = n_distinct(.taxonomy_key),
    n_abundance_profiles = n_distinct(.profile_key),
    .groups = "drop"
  )

conflicting_duplicates <- duplicate_check %>%
  filter(
    n_taxonomy_profiles > 1 |
      n_abundance_profiles > 1
  )

if (nrow(conflicting_duplicates) > 0) {
  readr::write_csv(
    conflicting_duplicates,
    file.path(output_dir, "ERROR_conflicting_duplicate_MAGs.csv")
  )
  stop(
    "同一 MAGs 存在不同分类或不同丰度记录，已停止分析；",
    "请查看 ERROR_conflicting_duplicate_MAGs.csv。"
  )
}

global_unique <- global_raw %>%
  distinct(MAGs, .keep_all = TRUE)

message(
  "全球丰度表：原始 ", nrow(global_raw),
  " 行；按 MAGs 去重后 ", nrow(global_unique), " 个唯一 MAG。"
)

# ---------- 7. 转成长表，并计算每个样品全部唯一 MAG 总丰度 ----------
global_long <- global_unique %>%
  select(
    MAGs, Family_clean, Genus_clean,
    all_of(sample_cols)
  ) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample_name",
    values_to = "Abundance"
  ) %>%
  mutate(Abundance = replace_na(as.numeric(Abundance), 0))

sample_total <- global_long %>%
  group_by(Sample_name) %>%
  summarise(
    Total_MAG_abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

if (any(sample_total$Total_MAG_abundance <= 0)) {
  bad_samples <- sample_total %>%
    filter(Total_MAG_abundance <= 0) %>%
    pull(Sample_name)
  stop(
    "以下样品的全部唯一 MAG 总丰度小于或等于 0：",
    paste(bad_samples, collapse = ", ")
  )
}

# ---------- 8. 分类层级分析函数 ----------
analyse_taxonomic_rank <- function(
    rank_col,
    local_rank_col,
    rank_label,
    file_prefix
) {

  local_taxon_map <- local_mag_keep %>%
    filter(!is.na(.data[[local_rank_col]])) %>%
    transmute(
      Taxon = .data[[local_rank_col]],
      MAG
    ) %>%
    group_by(Taxon) %>%
    summarise(
      Local_MAG_count = n_distinct(MAG),
      Local_MAG_IDs = collapse_unique(MAG),
      .groups = "drop"
    ) %>%
    arrange(Taxon)

  target_taxa <- local_taxon_map$Taxon

  global_taxa <- global_unique %>%
    pull(all_of(rank_col)) %>%
    unique() %>%
    na.omit() %>%
    as.character()

  # 严格名称匹配：只保留完全相同名称
  matched_taxa <- sort(intersect(target_taxa, global_taxa))
  unmatched_taxa <- sort(setdiff(target_taxa, global_taxa))

  message(
    rank_label, "：本地目标 ", length(target_taxa),
    " 个；严格匹配成功 ", length(matched_taxa),
    " 个；未匹配 ", length(unmatched_taxa), " 个。"
  )

  if (length(matched_taxa) == 0) {
    stop(rank_label, " 没有任何严格匹配成功的分类单元。")
  }

  # 构建“每个样品 × 每个匹配分类单元”的完整组合
  taxon_abundance_observed <- global_long %>%
    filter(.data[[rank_col]] %in% matched_taxa) %>%
    transmute(
      Sample_name,
      Taxon = .data[[rank_col]],
      Abundance
    ) %>%
    group_by(Sample_name, Taxon) %>%
    summarise(
      Taxon_abundance = sum(Abundance, na.rm = TRUE),
      .groups = "drop"
    )

  taxon_sample <- tidyr::expand_grid(
    Sample_name = sample_cols,
    Taxon = matched_taxa
  ) %>%
    left_join(
      taxon_abundance_observed,
      by = c("Sample_name", "Taxon")
    ) %>%
    mutate(
      Taxon_abundance = replace_na(Taxon_abundance, 0)
    ) %>%
    left_join(
      sample_total,
      by = "Sample_name",
      relationship = "many-to-one"
    ) %>%
    mutate(
      Taxon_relative_abundance_percent =
        Taxon_abundance / Total_MAG_abundance * 100
    )

  # 每个样品中所有匹配目标的累积相对丰度
  sample_summary <- taxon_sample %>%
    group_by(Sample_name) %>%
    summarise(
      Detected_target_taxa = sum(Taxon_abundance > 0),
      Detected_taxa_names = collapse_unique(
        Taxon[Taxon_abundance > 0]
      ),
      Target_abundance_sum = sum(Taxon_abundance),
      Total_MAG_abundance = first(Total_MAG_abundance),
      Target_relative_abundance_percent =
        Target_abundance_sum / Total_MAG_abundance * 100,
      Positive_sample = Target_abundance_sum > 0,
      .groups = "drop"
    ) %>%
    left_join(
      site_meta,
      by = "Sample_name",
      relationship = "one-to-one"
    ) %>%
    relocate(
      Sample_name, Site, Location_ID,
      longitude, latitude, sample_type,
      Region, Country
    )

  # 72 个 Site ID 的统计，仅作为核对，不直接决定地图点数
  site_summary <- sample_summary %>%
    group_by(
      Site, longitude, latitude,
      sample_type, Region, Country
    ) %>%
    summarise(
      sample_size = n_distinct(Sample_name),
      positive_sample_count = sum(Positive_sample),
      mean_relative_abundance_percent =
        mean(Target_relative_abundance_percent),
      min_relative_abundance_percent =
        min(Target_relative_abundance_percent),
      max_relative_abundance_percent =
        max(Target_relative_abundance_percent),
      Positive_site = any(Positive_sample),
      .groups = "drop"
    )

  # 地图使用“唯一经纬度”聚合。
  # 当前数据有 72 个 Site ID，但只有 52 对唯一经纬度。
  location_summary <- sample_summary %>%
    group_by(
      Location_ID, longitude, latitude,
      sample_type, Region, Country
    ) %>%
    summarise(
      Site_ID_count = n_distinct(Site),
      Site_IDs = collapse_unique(Site),
      sample_size = n_distinct(Sample_name),
      positive_sample_count = sum(Positive_sample),
      mean_relative_abundance_percent =
        mean(Target_relative_abundance_percent),
      min_relative_abundance_percent =
        min(Target_relative_abundance_percent),
      max_relative_abundance_percent =
        max(Target_relative_abundance_percent),
      Positive_location = any(Positive_sample),
      .groups = "drop"
    ) %>%
    arrange(longitude, latitude)

  # 每个分类单元在样品、Site ID、唯一地理点中的分布
  taxon_sample_meta <- taxon_sample %>%
    left_join(
      site_meta,
      by = "Sample_name",
      relationship = "many-to-one"
    )

  taxon_location <- taxon_sample_meta %>%
    group_by(
      Taxon, Location_ID,
      longitude, latitude
    ) %>%
    summarise(
      Mean_location_relative_abundance_percent =
        mean(Taxon_relative_abundance_percent),
      Positive_location = any(Taxon_abundance > 0),
      .groups = "drop"
    )

  taxon_summary <- taxon_sample_meta %>%
    group_by(Taxon) %>%
    summarise(
      Positive_samples = sum(Taxon_abundance > 0),
      Positive_Site_IDs =
        n_distinct(Site[Taxon_abundance > 0]),
      Positive_geographic_locations =
        n_distinct(Location_ID[Taxon_abundance > 0]),
      Min_positive_sample_relative_abundance_percent =
        positive_min(Taxon_relative_abundance_percent),
      Max_sample_relative_abundance_percent =
        positive_max(Taxon_relative_abundance_percent),
      .groups = "drop"
    ) %>%
    left_join(
      taxon_location %>%
        group_by(Taxon) %>%
        summarise(
          Min_positive_location_mean_percent =
            positive_min(
              Mean_location_relative_abundance_percent
            ),
          Max_location_mean_percent =
            positive_max(
              Mean_location_relative_abundance_percent
            ),
          .groups = "drop"
        ),
      by = "Taxon"
    ) %>%
    left_join(
      local_taxon_map,
      by = "Taxon",
      relationship = "one-to-one"
    ) %>%
    relocate(
      Taxon, Local_MAG_count, Local_MAG_IDs
    )

  unmatched_summary <- local_taxon_map %>%
    filter(Taxon %in% unmatched_taxa) %>%
    mutate(
      Positive_samples = 0L,
      Positive_Site_IDs = 0L,
      Positive_geographic_locations = 0L,
      Reason = paste0(
        "全球唯一 MAG 数据中不存在完全相同的 ",
        rank_label,
        " 名称；不使用其他分类层级替代。"
      )
    )

  # ---------- 输出统计表 ----------
  readr::write_csv(
    sample_summary,
    file.path(
      output_dir,
      paste0(file_prefix, "_sample_level_193_samples.csv")
    )
  )

  readr::write_csv(
    site_summary,
    file.path(
      output_dir,
      paste0(file_prefix, "_site_ID_level_72_sites.csv")
    )
  )

  readr::write_csv(
    location_summary,
    file.path(
      output_dir,
      paste0(file_prefix, "_geographic_location_level.csv")
    )
  )

  readr::write_csv(
    taxon_summary,
    file.path(
      output_dir,
      paste0(file_prefix, "_taxon_summary.csv")
    )
  )

  readr::write_csv(
    unmatched_summary,
    file.path(
      output_dir,
      paste0(file_prefix, "_strict_unmatched_taxa.csv")
    )
  )

  # ---------- 绘制参考文献风格全球地图 ----------
  # map_data() 是 ggplot2 的函数；maps 包只提供底层地图数据库。
  world_map <- ggplot2::map_data("world")

  plot_data <- location_summary %>%
    filter(Positive_location) %>%
    mutate(
      Carrier = "plastisphere"
    )

  readr::write_csv(
    plot_data %>%
      select(
        Location_ID, longitude, latitude,
        sample_type, Region, Country,
        Site_ID_count, Site_IDs, sample_size,
        mean_relative_abundance_percent
      ),
    file.path(
      output_dir,
      paste0(file_prefix, "_EXACT_coordinates_used_for_map.csv")
    )
  )

  base_plot <- ggplot() +
    geom_polygon(
      data = world_map,
      aes(x = long, y = lat, group = group),
      fill = "#EEEEEE",
      color = "#A6A6A6",
      linewidth = 0.22
    )

  if (bubble_size_by == "sample_count") {

    map_plot <- base_plot +
      geom_point(
        data = plot_data,
        aes(
          x = longitude,
          y = latitude,
          size = sample_size,
          fill = mean_relative_abundance_percent,
          shape = Carrier
        ),
        color = "#777777",
        stroke = 0.45,
        alpha = 0.88
      ) +
      scale_size_continuous(
        name = "No. of samples",
        range = c(2.2, 10),
        breaks = c(1, 3, 6, 9, 12),
        limits = c(1, 12)
      )

  } else if (bubble_size_by == "abundance") {

    map_plot <- base_plot +
      geom_point(
        data = plot_data,
        aes(
          x = longitude,
          y = latitude,
          size = mean_relative_abundance_percent,
          fill = mean_relative_abundance_percent,
          shape = Carrier
        ),
        color = "#777777",
        stroke = 0.45,
        alpha = 0.88
      ) +
      scale_size_continuous(
        name = paste0(
          rank_label,
          " relative abundance (%)"
        ),
        range = c(2.2, 10)
      )

  } else {
    stop(
      "bubble_size_by 只能设置为 'sample_count' 或 'abundance'。"
    )
  }

  map_plot <- map_plot +
    scale_shape_manual(
      name = "Carrier",
      values = c(plastisphere = 21)
    ) +
    scale_fill_gradientn(
      colours = bubble_palette,
      name = paste0(
        "Mean ", rank_label,
        "\nrelative abundance (%)"
      ),
      trans = scales::sqrt_trans(),
      labels = scales::label_number(
        accuracy = 0.001,
        trim = TRUE
      )
    ) +
    scale_x_continuous(
      breaks = seq(-180, 180, by = 60),
      labels = format_lon,
      limits = c(-180, 180),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq(-90, 90, by = 30),
      labels = format_lat,
      limits = c(-90, 90),
      expand = c(0, 0)
    ) +
    coord_fixed(
      ratio = 1,
      xlim = c(-180, 180),
      ylim = c(-90, 90),
      expand = FALSE,
      clip = "on"
    ) +
    labs(
      title = paste0(
        "Global distribution of MAG-associated ",
        tolower(rank_label),
        " in the plastisphere"
      ),
      subtitle = paste0(
        "Strict ", rank_label,
        " matching; one bubble per exact longitude–latitude pair"
      ),
      x = "Longitude",
      y = "Latitude",
      caption = paste0(
        nrow(plot_data),
        " geographic locations; ",
        sum(plot_data$sample_size),
        " samples represented"
      )
    ) +
    guides(
      shape = guide_legend(
        order = 1,
        override.aes = list(
          size = 3,
          fill = "white",
          alpha = 1
        )
      ),
      size = guide_legend(
        order = 2,
        override.aes = list(
          fill = "#999999",
          alpha = 1
        )
      ),
      fill = guide_colourbar(
        order = 3,
        barheight = grid::unit(32, "mm"),
        barwidth = grid::unit(5, "mm")
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.85
      ),
      panel.grid.major = element_line(
        color = "#C9CDD2",
        linewidth = 0.45
      ),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.55
      ),
      axis.ticks.length = grid::unit(2.5, "mm"),
      axis.text = element_text(
        color = "black",
        size = 10.5
      ),
      axis.title.x = element_text(
        face = "bold",
        size = 12,
        margin = margin(t = 8)
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 12,
        margin = margin(r = 8)
      ),
      plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 10,
        color = "#5F6368",
        hjust = 0.5
      ),
      plot.caption = element_text(
        size = 9,
        color = "#5F6368",
        hjust = 1
      ),
      legend.position = "left",
      legend.box = "vertical",
      legend.title = element_text(
        face = "bold",
        size = 10
      ),
      legend.text = element_text(
        size = 9
      ),
      plot.margin = margin(
        t = 10, r = 15, b = 8, l = 10
      )
    )

  # PNG
  ggsave(
    filename = file.path(
      output_dir,
      paste0(file_prefix, "_global_map.png")
    ),
    plot = map_plot,
    width = 14,
    height = 7.8,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  # 可供 Illustrator 编辑的矢量 PDF
  ggsave(
    filename = file.path(
      output_dir,
      paste0(file_prefix, "_global_map.pdf")
    ),
    plot = map_plot,
    width = 14,
    height = 7.8,
    units = "in",
    device = "pdf",
    useDingbats = FALSE,
    bg = "white"
  )

  # RDS 便于后续直接修改图形对象
  saveRDS(
    map_plot,
    file.path(
      output_dir,
      paste0(file_prefix, "_global_map_plot.rds")
    )
  )

  list(
    rank_label = rank_label,
    target_taxa = target_taxa,
    matched_taxa = matched_taxa,
    unmatched_taxa = unmatched_taxa,
    sample_summary = sample_summary,
    site_summary = site_summary,
    location_summary = location_summary,
    taxon_summary = taxon_summary,
    unmatched_summary = unmatched_summary,
    plot = map_plot
  )
}

# ---------- 9. 严格 Genus 匹配 ----------
genus_result <- analyse_taxonomic_rank(
  rank_col = "Genus_clean",
  local_rank_col = "Genus",
  rank_label = "Genus",
  file_prefix = "strict_genus"
)

# ---------- 10. 严格 Family 匹配 ----------
family_result <- analyse_taxonomic_rank(
  rank_col = "Family_clean",
  local_rank_col = "Family",
  rank_label = "Family",
  file_prefix = "strict_family"
)

# ---------- 11. 输出总体核对信息 ----------
overall_summary <- tibble(
  Analysis = c("Strict Genus", "Strict Family"),
  Local_target_taxa = c(
    length(genus_result$target_taxa),
    length(family_result$target_taxa)
  ),
  Strict_matched_taxa = c(
    length(genus_result$matched_taxa),
    length(family_result$matched_taxa)
  ),
  Strict_unmatched_taxa = c(
    length(genus_result$unmatched_taxa),
    length(family_result$unmatched_taxa)
  ),
  Positive_samples = c(
    sum(genus_result$sample_summary$Positive_sample),
    sum(family_result$sample_summary$Positive_sample)
  ),
  Total_samples = n_distinct(site_meta$Sample_name),
  Positive_Site_IDs = c(
    sum(genus_result$site_summary$Positive_site),
    sum(family_result$site_summary$Positive_site)
  ),
  Total_Site_IDs = n_distinct(site_meta$Site),
  Positive_geographic_locations = c(
    sum(genus_result$location_summary$Positive_location),
    sum(family_result$location_summary$Positive_location)
  ),
  Total_geographic_locations =
    n_distinct(site_meta$Location_ID),
  Min_positive_sample_relative_abundance_percent = c(
    positive_min(
      genus_result$sample_summary$
        Target_relative_abundance_percent
    ),
    positive_min(
      family_result$sample_summary$
        Target_relative_abundance_percent
    )
  ),
  Max_sample_relative_abundance_percent = c(
    positive_max(
      genus_result$sample_summary$
        Target_relative_abundance_percent
    ),
    positive_max(
      family_result$sample_summary$
        Target_relative_abundance_percent
    )
  ),
  Min_positive_location_mean_percent = c(
    positive_min(
      genus_result$location_summary$
        mean_relative_abundance_percent
    ),
    positive_min(
      family_result$location_summary$
        mean_relative_abundance_percent
    )
  ),
  Max_location_mean_percent = c(
    positive_max(
      genus_result$location_summary$
        mean_relative_abundance_percent
    ),
    positive_max(
      family_result$location_summary$
        mean_relative_abundance_percent
    )
  )
)

readr::write_csv(
  overall_summary,
  file.path(output_dir, "overall_genus_family_summary.csv")
)

message("分析完成。输出目录：", normalizePath(output_dir))
message(
  "经纬度核验：",
  n_distinct(site_meta$Sample_name), " 个样品；",
  n_distinct(site_meta$Site), " 个 Site ID；",
  n_distinct(site_meta$Location_ID), " 对唯一经纬度。"
)
