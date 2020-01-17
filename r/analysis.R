################################################################################
## Code for manuscript: "Catalyzing sustainable fisheries management though behavior change interventions"
## Authors: Gavin McDonald, Molly Wilson, Diogo Veríssimo, Rebecca Twohey, Michaela Clemence, Dean Apistar, Stephen Box, Paul Butler, Fel Cesar Cadiz, Stuart J. Campbell, Courtney Cox, Micah Effron, Steve Gaines, Raymond Jakub, Roquelito H. Mancao, Pablo T. Rojas, Rocky Sanchez Tirona, Gabriel Vianna
## Contact: Gavin McDonald (gmcdonald@bren.ucsb.edu)
################################################################################

################################################################################
## Setup
################################################################################

# Load packages
library(tidyverse)
library(stringr)
library(cowplot)
library(broom)
library(sandwich)
library(lmtest)
library(lavaan)
library(semPlot)
library(janitor)
library(rworldmap)
library(rmapshaper)
library(sf)
library(corrplot)
library(here)
# Load data
monitoring_data <- read_csv(here::here("data/monitoring_data.csv")) %>%
  # Don't include habitat data for paper
  filter(survey != "Sustainable Habitat")

matching_data <- read_csv("data/phils_matching_data.csv") %>%
  group_by(site,control_impact,country,survey) %>%
  summarize(matching_score = sum(value)) %>%
  ungroup() %>%
  mutate(control_impact = as.numeric(control_impact))

site_data <- read_csv(here::here("data/site_data.csv"))

survey_questions <- read_csv(here::here("data/survey_question_lookup.csv"))


# Proess data
# Define which condensed indicators belong in each survey
kap_indicators <- c("Knowledge","Attitude","Communication")
bc_indicators <- c("Catch Reporting","Enforcement","Licensing","TURF Compliance","NTZ Compliance","Management Participation")
eco_indicators <- c("Biomass Inside NTZ","Biomass Outside NTZ")
social_indicators <- c("Food Security","Livelihood Stability","Political Trust",
                       "Social Equity","Reduced Fishing Costs","Social Trust",
                       "Subjective Well-being","Catch Trend 5 Years",
                       "Catch Trend 1 Year","Collective Efficacy","Household Assets","Reduced Fishing Costs")

# For the social data, create composite scores for community support, sustainable fishing practices, and sustainable livelihoods indicators
# Sometimes, multiple questions may be asked relating to a single category to each person
# Take average of similar questions to create composite score
monitoring_data_composite_full <- monitoring_data %>%
  filter(survey %in% c("Community Support","Sustainable Fishing Practices","Sustainable Livelihoods")) %>%
  mutate(indicator = indicator_full) %>%
  group_by(country,site,survey,before_after,control_impact,replicate,indicator) %>%
  summarize(value = mean(value,na.rm=TRUE)) %>%
  ungroup() %>%
  # Add back in other indicators
  bind_rows(monitoring_data %>%
    filter(!survey %in% c("Community Support","Sustainable Fishing Practices","Sustainable Livelihoods")) %>%
      mutate(indicator = indicator_full))

# Do the same, but using condensed indicators for knowledge, attitude, and communication
monitoring_data_composite_condensed <- monitoring_data %>%
  filter(survey %in% c("Community Support","Sustainable Fishing Practices","Sustainable Livelihoods")) %>%
  mutate(indicator = indicator_condensed) %>%
  group_by(country,site,survey,before_after,control_impact,replicate,indicator) %>%
  summarize(value = mean(value,na.rm=TRUE)) %>%
  ungroup() %>%
  # Add back in other indicators
  bind_rows(monitoring_data %>%
              filter(!survey %in% c("Community Support","Sustainable Fishing Practices","Sustainable Livelihoods")) %>%
              mutate(indicator = indicator_condensed))

################################################################################
## Code for Before-after and diff-in-diff Analysis analysis
## Code written by Gavin McDonald (gmcdonald@brebn.ucsb.edu)
## November 2, 2018
################################################################################

# Custom function for linear model that is robust to heteroskedasticity
# Scale value to get standardized beta coefficients
lm_robust <- function(formula,df){
  lm(formula, data = df %>%
       mutate(value = scale(value))) %>%
    lmtest::coeftest(vcov = sandwich::vcovHC(., "HC1"))
}

# Before-after analysis for site-indicator pairs
model_beforeafter_site <- monitoring_data_composite_condensed %>%
  # Only look at FF sites
  filter(control_impact == 1) %>%
  group_by(country,site,indicator,survey) %>%
  # Use custom lm model to include heteroskedasticity
  do(before_after_model = safely(lm_robust)("value ~ before_after", df = .)$result) %>% 
  # # Use 95% confidence intervals
  tidy(before_after_model) %>%
  # tidy(did_model,conf.int = TRUE, conf.level = 0.95) %>%
  mutate(conf.high.95 = estimate + 1.96 * std.error,
         conf.low.95 = estimate - 1.96 * std.error) %>%
  # Only take standardized before after estimaor
  filter(term == "before_after") %>%
  ungroup()

# Diff-in-diff analysis for site-indicator pairs
model_DiD_site<- monitoring_data_composite_condensed %>%
  # Add matching data
  left_join(matching_data, by = c("site","control_impact","country","survey")) %>%
  group_by(country, site,indicator,survey) %>%
  # Saturated model. Include interaction to look at causality
  # Include matching score
  # Use custom lm model to include heteroskedasticity
  do(before_after_model = safely(lm_robust)("value ~ before_after * control_impact + matching_score", df = .)$result) %>% 
  # # Use 95% confidence intervals
  tidy(before_after_model,conf.int = TRUE, conf.level = 0.95) %>%
  # tidy(did_model,conf.int = TRUE, conf.level = 0.95) %>%
  mutate(conf.high.95 = estimate + 1.96 * std.error,
         conf.low.95 = estimate - 1.96 * std.error) %>%
  # # Ignore the intercept estimates
  filter(term == "before_after:control_impact") %>%
  ungroup()

# Join before-after and diff-in-diff models
site_models <- bind_rows(model_beforeafter_site,model_DiD_site)

# Calculate precision-weight estimates for each country-indicator pair across sites
# Follow this procedure: https://www.ucdmc.ucdavis.edu/ctsc/area/Resource_Library/documents/Stewart%20-%20Meta-analysis%20-2017.png
country_precision_weighted_estimates <- site_models %>%
  # Remove empty values, or those with 0 standard error
  filter(!is.nan(estimate) & !is.nan(std.error) & std.error != 0) %>%
  #mutate(country = ifelse(term == "before_after","Aggregate",country)) %>%
  mutate(weight = 1 / std.error ^2,
         weighted_estimate = weight * estimate) %>%
  group_by(country,survey,indicator,term) %>%
  summarize(weighted_estimate = sum(weighted_estimate) / sum(weight),
            weight_estimate_standard_error = 1 / sqrt(sum(weight))) %>%
  mutate(conf.low.95 = weighted_estimate - 1.96 * weight_estimate_standard_error,
         conf.high.95 = weighted_estimate + 1.96 * weight_estimate_standard_error) %>%
  ungroup()

################################################################################
## Figures
################################################################################

# Set figure color palette


palPositive<- setNames(c("forestgreen","darkorchid","white",
                         "forestgreen","darkorchid","white"),
                       c("Positive\n(p-value < 0.05)",
                         "Negative\n(p-value < 0.05)",
                         "Not Significant\n(p-value >= 0.05)",
                         "Positive\n(95% C.I. range is entirely positive)",
                         "Negative\n(95% C.I. range is entirely negative)",
                         "Not Significant\n(95% C.I. range overlaps 0)"))
palSignif <- setNames(c(0.1,1),c("No","Yes"))

# Custom function to make paper figures
paper_figure <- function(indicator_filter, term_filter, file_name,fig_height,site_models_fig,pw_models_fig){
  
  # Filtered precision-weighted estimates dataframe
  temp_df_weighted <- pw_models_fig %>%
    filter(!! indicator_filter) %>%
    filter(!! term_filter) %>%
    mutate(survey = case_when(survey == "Community Support" ~ "Community\nSupport",
                              survey ==  "Sustainable Fishing Practices" ~ "Sustainable Fishing\nPractices",
                              survey == "Sustainable Livelihoods" ~ "Sustainable Livelihoods",
                              survey ==  "Sustainable Ecosystem" ~ "Sustainable\nEcosystem"))
  
  # Filtered site dataframe
  temp_df_site <- site_models_fig %>%
    filter(!! indicator_filter) %>%
    filter(!! term_filter)  %>%
    mutate(survey = case_when(survey == "Community Support" ~ "Community\nSupport",
                              survey ==  "Sustainable Fishing Practices" ~ "Sustainable Fishing\nPractices",
                              survey == "Sustainable Livelihoods" ~ "Sustainable Livelihoods",
                              survey ==  "Sustainable Ecosystem" ~ "Sustainable\nEcosystem"))
  
  # Reorder facets for plotting
  temp_df_weighted$survey <- factor(temp_df_weighted$survey, 
                                    levels = c("Community\nSupport",
                                               "Sustainable Fishing\nPractices",
                                               "Sustainable\nEcosystem",
                                               "Sustainable Livelihoods"))
  
  temp_df_site$survey <- factor(temp_df_site$survey, 
                                levels = c("Community\nSupport",
                                           "Sustainable Fishing\nPractices",
                                           "Sustainable\nEcosystem",
                                           "Sustainable Livelihoods"))
  
  
  # Set figure limit the same for each plot
  x_range <- c(min(c(temp_df_site$estimate,temp_df_weighted$conf.low.95)),
               max(c(temp_df_site$estimate,temp_df_weighted$conf.high.95))) %>%
    extendrange(f = 0.05)
  
  # Panel B
  
  p1 <- temp_df_weighted %>%
    mutate(Change = case_when(conf.low.95 < 0 & conf.high.95 > 0 ~ "Not Significant\n(95% C.I. range overlaps 0)",
                              weighted_estimate > 0 ~ "Positive\n(95% C.I. range is entirely positive)",
                              TRUE ~ "Negative\n(95% C.I. range is entirely negative)")) %>%
    ggplot(aes(x = indicator, y = weighted_estimate, fill = `Change`)) +
    geom_hline(yintercept = 0) +
    geom_errorbar(aes(ymin = conf.low.95,
                      ymax = conf.high.95),
                  width = 0.5) +
    geom_point(size = 3,shape=21,color="black") +
    coord_flip() +
    xlab("") +
    ylim(x_range) +
    ylab("Effect Size") +
    scale_fill_manual(values = palPositive) +
    theme_bw(base_size = 12) +
    facet_grid(survey~country,scales="free", space="free_y",switch="y") +
    theme(axis.text.y = element_blank(),
          strip.background.y = element_blank(),
          strip.text.y = element_blank(),
          strip.background=element_rect(fill="white"),
          legend.direction = "vertical", 
          legend.position = "bottom",
          legend.box = "horizontal") + 
    guides(alpha = guide_legend(override.aes = list(fill = "black")))
    #guides(colour=FALSE,alpha=FALSE,fill=FALSE)
  
  indicator_order <- layer_scales(p1)$x$range$range
  
  # Panel A
  
  p2 <- temp_df_site %>%
    mutate(Change = case_when(p.value >= 0.05 | is.na(p.value) ~ "Not Significant\n(p-value >= 0.05)",
                              estimate > 0 ~ "Positive\n(p-value < 0.05)",
                              TRUE ~ "Negative\n(p-value < 0.05)")) %>%
    ggplot(aes(x = estimate, y = indicator)) +
    geom_vline(xintercept = 0) +
    geom_jitter(height = 0.2,aes(fill = `Change`), colour="black",shape =21,size=3) +
    ylab("") +
    xlim(x_range) +
    xlab("Effect Size") +
    scale_fill_manual(values = palPositive) +
    theme_bw(base_size = 12) +
    facet_grid(survey~country,scales="free", space="free_y",switch="y") +
    theme(legend.direction = "vertical", 
          legend.position = "bottom",
          legend.box = "horizontal",
          strip.background=element_rect(fill="white")) + 
    guides(alpha = guide_legend(override.aes = list(fill = "black")))
  
  # Use cowplot to setup A and B panels
  p <- plot_grid(plot_grid(p2, p1, align = "h",rel_widths = c(1.75,1), labels = c("A", "B")),ncol=1, rel_heights = c(1, .1))
  # Save figure
  p
  ggsave(here::here(paste0("output_figures/",file_name,".png")),p,width = 7.5,height=fig_height,device="png",dpi=300)
}

# Make before-after figure  - All surveys
paper_figure(indicator_filter = expr(TRUE),
             term_filter = expr(term == "before_after" ),
             file_name = "figure_3_before_after",
             fig_height = 10,
             site_models_fig = site_models,
             pw_models_fig = country_precision_weighted_estimates)

# Diff-in-diff for main body of paper - Phils Sustainable Ecosystem and Sustainable Livelihoods
paper_figure(indicator_filter = expr((survey == "Sustainable Livelihoods" & country == "Philippines") |
                                       (survey == "Sustainable Ecosystem" & country == "Philippines")),
             term_filter = expr(term == "before_after:control_impact" ),
             file_name = "figure_4_diff_in_diff",
             fig_height = 6,
             site_models_fig = site_models,
             pw_models_fig = country_precision_weighted_estimates)

# Site summary figure
site_summary_fig <- site_data %>%
  mutate(percentage_reserve = area_current_reserves_ha / area_current_turfs_ha * 100) %>%
  filter(site_type == "Intervention") %>%
  rename(`Number of fishers` = fishers_ff_target_communities,
                `TURF area [Ha]` = area_current_turfs_ha,
                `NTZ area [Ha]` = area_current_reserves_ha,
                `Percentage TURF covered by NTZ [%]` = percentage_reserve) %>%
  filter(`NTZ area [Ha]`<`TURF area [Ha]`) %>%
  mutate(`Number of fishers` = as.numeric(`Number of fishers`)) %>%
  gather(indicator,value,`Number of fishers`:`Percentage TURF covered by NTZ [%]` ) %>%
  ggplot() +
  geom_boxplot(aes(x=country,y=value),outlier.alpha = 0)+
  geom_jitter(aes(x=country,y=value),width=0.1,shape=21,size=1) +
  facet_wrap(.~indicator,scales="free",strip.position = "left")+
  theme_bw(base_size=12) + 
  scale_y_continuous(labels=scales::comma,trans="log")+
  labs(x="",
       y = "") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill=NA))

ggsave(here::here("output_figures/figure_s1_site_summary.png"),site_summary_fig,width = 7.5,height=7.5,device="png",dpi=300)

################################################################################
## Tables
################################################################################

# Custom function for matching table
matching_table_func <- function(indicator_filter,country,survey,matching_data_table){
  matching_data_table %>% 
    filter(!! indicator_filter) %>%
    left_join(site_data %>%
                dplyr::select(site_name,site), by = "site") %>%
    dplyr::select(-site) %>%
    rename(Country = country,
           Site = site_name,
           Survey = survey) %>%
    mutate(Survey = survey) %>%
    spread(control_impact,matching_score) %>%
    rename(`Control Site Score` = `0`,
           `Intervention Site Score` = `1`)%>%
    bind_rows(
      summarise_at(.,c("Control Site Score",
                       "Intervention Site Score"), mean, na.rm = TRUE) %>%
        mutate(Site = "Mean",
               Country = country,
               Survey = survey)
    ) %>%
    mutate(`Intervention Site Score` = signif(`Intervention Site Score`,3),
           `Control Site Score` = signif(`Control Site Score`,3)) %>%
    dplyr::select(Site,`Intervention Site Score`,`Control Site Score`)
}

# Diff-in-diff matching table for main-body - Phils sustainable ecosystem and livelihoods
matching_table_func(indicator_filter = expr(country == "Philippines" &
                                              survey %in% c("Sustainable Livelihoods") & site %in% c("CULASI_","DAPACOR","GUBATRA")),
                    country = "Philippines",
                    survey = "Sustainable Livelihoods",
                    matching_data_table = matching_data) %>%
  write_csv(path=here::here("output_tables/matching_table.csv"))



survey_order <- c("Community Support","Sustainable Fishing Practices","Sustainable Ecosystem","Sustainable Livelihoods")
survey_order_short <- c("CS","SFP","SE","SL")
site_models %>%
  mutate(term = case_when(term == "before_after" ~ "n Intervention",
                          term == "before_after:control_impact" ~ "n Control")) %>%
  group_by(survey,indicator,term,country) %>%
  summarize(count = n())  %>%
  ungroup() %>%
  gather(variable, value, -c(term,country,survey,indicator)) %>%
  unite(temp, country, term, sep=" ") %>%
  spread(temp, value) %>%
  dplyr::select(-variable) %>%
  replace(is.na(.), 0) %>%
  rename(Survey = survey,Indicator = indicator) %>%
  filter(!is.na(Survey))%>%
  mutate(Survey =  factor(Survey, levels = survey_order)) %>%
  arrange(Survey) %>%
  mutate(Survey = case_when(Indicator %in% kap_indicators ~ "CS",
                            Indicator %in% bc_indicators ~ "SFP",
                            Indicator %in% social_indicators ~"SL",
                            Indicator %in% eco_indicators ~ "SE"))%>%
  write_csv(path=here::here("output_tables/data_summary.csv"))

# Make tables that have sample sites for all sites and indicators

ss_table <- monitoring_data_composite_condensed %>% 
  group_by(country,site,indicator,before_after,control_impact) %>% 
  summarize(sample_size = n_distinct(replicate)) %>% 
  group_by(country,site,indicator) %>% 
  arrange(control_impact,before_after) %>%
  summarize(sample_size = paste(sample_size,collapse=", ")) %>% 
  ungroup() %>% 
  mutate(sample_size = sub("^((?:.*?,){1}.*?),", "\\1;", sample_size, perl = TRUE))%>%
  rename(Country = country,
         Site = site,
         Indicator = indicator) %>%
  mutate(Survey = case_when(Indicator %in% kap_indicators ~ "CS",
                            Indicator %in% bc_indicators ~ "SFP",
                            Indicator %in% social_indicators ~"SL",
                            Indicator %in% eco_indicators ~ "SE"))


ss_table_brazil <- ss_table %>%
  filter(Country == "Brazil")%>%
  spread(Site,sample_size)%>%
  mutate(Survey =  factor(Survey, levels = survey_order_short)) %>%
  arrange(Survey) %>%
  dplyr::select(Survey,Indicator,everything()) %>%
  dplyr::select(-Country)

ss_table_indo <- ss_table %>%
  filter(Country == "Indonesia")%>%
  spread(Site,sample_size)%>%
  mutate(Survey =  factor(Survey, levels = survey_order_short)) %>%
  arrange(Survey) %>%
  dplyr::select(Survey,Indicator,everything()) %>%
  dplyr::select(-Country)

ss_table_phils <- ss_table %>%
  filter(Country == "Philippines")%>%
  spread(Site,sample_size)%>%
  mutate(Survey =  factor(Survey, levels = survey_order_short)) %>%
  arrange(Survey) %>%
  dplyr::select(Survey,Indicator,everything()) %>%
  dplyr::select(-Country)

# Need to break tables up to fit into paper
write_csv(ss_table_brazil,path=here::here("output_tables/ss_table_brazil.csv"))
write_csv(ss_table_indo[,1:9],path=here::here("output_tables/ss_table_indo_1.csv"))
write_csv(ss_table_indo[,c(1,2,seq(10,16))],path=here::here("output_tables/ss_table_indo_2.csv"))
write_csv(ss_table_phils[,1:9],path=here::here("output_tables/ss_table_phils_1.csv"))
write_csv(ss_table_phils[,c(1,2,seq(10,16))],path=here::here("output_tables/ss_table_phils_2.csv"))
write_csv(ss_table_phils[,c(1,2,seq(17,23))],path=here::here("output_tables/ss_table_phils_3.csv"))

# Survey questions

survey_question_table_full <- survey_questions %>%
  dplyr::select(Country = country,
                Survey = survey,
                Indicator = indicator_full,
                Question = full_question,
                Responses = response_options) %>%
  distinct() %>%
  mutate(Question = str_to_title(Question)) %>%
  arrange(Country,Survey,Indicator)

survey_table_combinations <- survey_question_table_full %>%
  dplyr::select(Country,Survey) %>%
  distinct() %>%
  mutate(Country = as.list(Country),
         Survey = as.list(Survey))

make_table <- function(country,survey,df){
  tmp_table <- df %>%
    filter(Country == country,
           Survey == survey)
  write_csv(tmp_table,path=here::here(paste0("output_tables/","survey_questions_",country,"_",survey,".csv")) %>%
              tolower() %>%
              str_replace(" ","_"))
}

map2(survey_table_combinations$Country,
     survey_table_combinations$Survey,
     ~make_table(.x,.y,
                 survey_question_table_full))

# Table describing which sites obtained single indicators through multiple questions
multiple_question_table <- survey_questions %>%
  mutate(question = paste0("[Question: ",full_question,"; Response: ",response_options)) %>%
  group_by(country, site_code,survey,indicator_condensed) %>%
  mutate(question = paste0(row_number(),": ",question)) %>%
  summarize(number_questions = n(),
            questions = paste(question,collapse="; ")) %>%
  filter(number_questions > 1)

write_csv(multiple_question_table,path=here::here("output_tables/multiple_question_table.csv"))

################################################################################
## Code for SEM analysis in manuscript
## Code written by Molly Wilson (mwwilson@ucsb.edu)
## November 2, 2018
################################################################################

# Convert indicators into columns 
dat_sem <- monitoring_data_composite_condensed %>%
  # Select only impact sites
  filter(control_impact == 1) %>% 
  # Calculate mean before and after values for each indicator at each site
  group_by(country,site,indicator,before_after) %>%
  summarize(mean_value = mean(value,na.rm=TRUE)) %>%
  ungroup() %>%
  # Calculate difference in before/after means for each indicator at each site
  spread(before_after,mean_value) %>%
  mutate(before_after_diff = `1` - `0`) %>%
  filter(!is.nan(before_after_diff) & !is.na(before_after_diff)) %>%
  dplyr::select(-`0`,-`1`) %>%
  spread(indicator,before_after_diff) %>% 
  clean_names(case = "old_janitor") %>%
  rename(catch_rep=catch_reporting, 
         mgmt_part=management_participation, 
         ntz_comp=ntz_compliance, 
         turf_comp=turf_compliance, 
         bm_in=biomass_inside_ntz, 
         bm_out=biomass_outside_ntz, 
         catch_tr1=catch_trend_1_year, 
         catch_tr5=catch_trend_5_years,  
         col_eff=collective_efficacy,
         assets=household_assets, 
         live_stab=livelihood_stability,
         pol_trust=political_trust, 
         sub_wellbeing=subjective_well_being)


################################################################################
## SEM Analysis
################################################################################

# Running two SEM model iterations based on Theory of Change relationships, with and without the inclusion of ecological variables

# Excluding ecological variables:
# Define model pathways
toc1 <- '
      # latent variable indicators  
          com_support =~ attitude + knowledge + communication
      # regressions
          mgmt_part ~ com_support
          enforcement ~ com_support
          ntz_comp ~ com_support + mgmt_part + enforcement
          turf_comp ~ com_support + mgmt_part + enforcement
          catch_rep ~ com_support
    '

# Run and inspect SEM:
fit1 <- sem(toc1, data=dat_sem, missing="ML", std.lv=TRUE)
summary(fit1, fit.measures=TRUE, standardized=TRUE)
fitMeasures(fit1, c("cfi","rmsea","srmr"))


# Including ecological variables:
# Define model pathways
toc2 <- '
      # latent variable indicators  
          com_support =~ attitude + knowledge + communication
      # regressions
          mgmt_part ~ com_support
          enforcement ~ com_support
          ntz_comp ~ com_support + mgmt_part + enforcement
          turf_comp ~ com_support + mgmt_part + enforcement
          catch_rep ~ com_support
          bm_in ~ ntz_comp
          bm_out ~ turf_comp + bm_in
    '

# run and inspect SEM:
fit2 <- sem(toc2, data=dat_sem, missing="ML", std.lv=TRUE)
summary(fit2, fit.measures=TRUE, standardized=T)
fitMeasures(fit2, c("cfi","rmsea","srmr"))


################################################################################
## SEM Figures
################################################################################

sem1 <- semPaths(fit1, 'std', curvePivot=TRUE, layout = "tree", intercept=FALSE, residuals = FALSE, nCharNodes = 0, 
                 label.scale=FALSE, sizeMan = 8, sizeLat = 10, label.cex=.4, esize=.5)
sem2 <- semPaths(fit2, 'std', curvePivot=TRUE, layout = "tree2", intercept=FALSE, residuals = FALSE, nCharNodes = 0, 
                 label.scale=FALSE, sizeMan = 8, sizeLat = 10, label.cex=.4, esize=.5)

################################################################################
## Map Figure
################################################################################

# Load world land map, use robinson projection
world_land <- st_as_sf(rworldmap::countriesLow) %>%
  st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# Creat world bounding box
world_bbox <- st_polygon(list(rbind(c(-88.5,-38.9),c(-88.5,28.9),c(156.9,28.9),c(156.9,-38.9),c(-88.5,-38.9)))) %>% st_sfc()%>% 
  `st_crs<-`("+proj=longlat +datum=WGS84 +no_defs") %>%
  st_transform(st_crs(world_land)) %>% st_bbox()

site_coords <- site_data %>%
  filter(!is.na(lat) & !is.na(lon)) %>%
  st_as_sf(coords = c("lon","lat")) %>% 
  `st_crs<-`("+proj=longlat +datum=WGS84 +no_defs") %>%
  st_transform(st_crs(world_land))

site_coords$site_type <- factor(site_coords$site_type, levels = c("Intervention","Control"))



# Get bounding box corrdinate for each country by buffering 5 degrees around country outline
brazil_bbox <- world_land %>%
  filter(ISO3 == "BRA") %>%
  st_buffer(1) %>%
  st_bbox()

indonesia_bbox <- world_land %>%
  filter(ISO3 == "IDN") %>%
  st_buffer(1) %>%
  st_bbox()

philippines_bbox <- world_land %>%
  filter(ISO3 == "PHL") %>%
  st_buffer(1) %>%
  st_bbox()
# Map color palette
pal <- setNames(c("dodgerblue","gold3"),c("Intervention","Control"))
# Make shape palette
shape_pal <- setNames(c(18,17),c("Intervention","Control"))

# Build base map, which we'll use for rest of maps
base_map <- ggplot()+ 
  geom_sf(data=world_land %>%
            mutate(Country = ifelse(ISO3 %in% c("BRA","PHL","IDN"),"Target","Other")), 
          color = NA, aes(fill=Country)) +
  scale_fill_manual(values=c("grey85","grey60")) +
  scale_color_manual(values = pal,name = "Site Type") +
  scale_shape_manual(values = shape_pal,name = "Site Type") +
  geom_sf(data = site_coords, aes(shape=site_type),color = "black",size=2.5, show.legend = FALSE) +
  geom_sf(data = site_coords,aes(color=site_type,shape=site_type), show.legend = "point",size=2) +
  guides(fill=FALSE)+
  theme_bw() +
  theme(plot.title = element_text(color="black",hjust=0,vjust=1, size=rel(1)),
        panel.border = element_rect(color = 'black', fill = NA,size=1),
        panel.grid.major = element_line(colour = "white"), 
        legend.text = element_text(color = "black", size = rel(2)),
        legend.title = element_text(color = "black", size = rel(2)),
        legend.title.align = 1,
        legend.background = element_rect(fill="white"),
        legend.position = "right",
        legend.direction="vertical",
        legend.text.align = 0,
        axis.text = element_text(color = "black", size = rel(1)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  guides(colour = guide_legend(override.aes = list(size=4)))

world_map <- base_map +
  geom_sf(data = philippines_bbox%>%
            st_as_sfc(),fill=NA) +
  geom_sf(data = indonesia_bbox%>%
            st_as_sfc(),fill=NA) +
  geom_sf(data = brazil_bbox%>%
            st_as_sfc(),fill=NA) +
  guides(color = FALSE,shape=FALSE)+
  xlim(c(world_bbox[1], world_bbox[3])) +
  ylim(c(world_bbox[2], world_bbox[4]))

brazil_map <- base_map +
  xlim(c(brazil_bbox[1], brazil_bbox[3])) +
  ylim(c(brazil_bbox[2], brazil_bbox[4])) +
  guides(color = FALSE,shape=FALSE)

indonesia_map <- base_map +
  xlim(c(indonesia_bbox[1], indonesia_bbox[3])) +
  ylim(c(indonesia_bbox[2], indonesia_bbox[4])) +
  guides(color = FALSE,shape=FALSE)

philippines_map <- base_map +
  xlim(c(philippines_bbox[1], philippines_bbox[3])) +
  ylim(c(philippines_bbox[2], philippines_bbox[4])) +
  guides(color = FALSE,shape=FALSE)

# Put all maps together
combined_map <- plot_grid(plot_grid(world_map,labels = "A"),
                          plot_grid(brazil_map, philippines_map, labels = c("B", "C"), rel_widths = c(0.65,0.35)),
                          plot_grid(indonesia_map, get_legend(base_map), labels = c("D", ""), rel_widths = c(0.65,0.35)),
                          ncol=1,
                          rel_heights = c(1,0.8,0.5))

ggsave(paste0(here::here("output_figures/figure_1_map.png")),combined_map,width = 7.5,height=7.5,device="png",dpi=300)