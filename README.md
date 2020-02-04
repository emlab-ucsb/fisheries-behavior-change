# fisheries-behavior-change
This repository contains the code and data for reproducing McDonald, *et al.* 2019 (Conservation Biology): "Catalyzing sustainable fisheries management though behavior change interventions" (https://doi.org/10.1111/cobi.13475).

> **Title**: Catalyzing sustainable fisheries management though behavior change interventions

> **Authors**: Gavin McDonald, Molly Wilson, Diogo Veríssimo, Rebecca Twohey, Michaela Clemence, Dean Apistar, Stephen Box, Paul Butler, Fel Cesar Cadiz, Stuart J. Campbell, Courtney Cox, Micah Effron, Steve Gaines, Raymond Jakub, Roquelito H. Mancao, Pablo T. Rojas, Rocky Sanchez Tirona, Gabriel Vianna


> **Abstract:** Small-scale fisheries are an important livelihood and primary protein source for coastal communities in many of the poorest regions in the world, yet many suffer from overfishing, requiring effective and scalable management solutions. Positive ecological and socioeconomic responses to management typically lag behind immediate costs borne by fishers from fishing pressure reductions necessary for fisheries recovery. These short-term costs challenge the long-term success of these interventions. However, social marketing may increase perceptions of management benefits before ecological and socioeconomic benefits are fully realized, driving new social norms and ultimately long-term sustainable behavior change. Using ecological surveys and community-perceived measures of management support and socioeconomic conditions, we assess the impact of a standardized small-scale fisheries management intervention that was implemented across 41 sites in Brazil, Indonesia, and the Philippines. The intervention combines TURF-reserves (community-based Territorial Use Rights for Fishing coupled with no-take marine reserves) with locally-tailored social marketing behavior change campaigns. Leveraging data across diverse indicators, our results suggest that communities were developing new social norms and fishing more sustainably, even before long- term ecological and socioeconomic benefits of fisheries management had materialized. 

## Repository structure  

```
fisheries-behavior-change 
  |__ data
  |__ output_figures
  |__ r
  |__ output_tables
```

## Software

This analysis was performed in R. The script for fully reproducing the paper analysis can be found in `r/analysis.R`.

## Data

The `data` folder contains four input data files, with full metadata given below:  

* `monitoring_data.csv`: All monitoring data for sites from the Community Support, Sustainable Fishing Practices, Sustainable Ecosystems, or Sustainable Livelihoods impact surveys   
* `site_data.csv`: Descriptive demographics and statistics for each site  
* `survey_question_lookup.csv`: A lookup table that matches specific socioeconomic survey questions from each site to specific indicators     
* `phils_matching_data.csv`: Site attribute scores used for control site selection for the 3 sets of matched impact / control sites for the Philippines sustainable ecosystems and sustainable livelihoods survey  

### `monitoring_data.csv`

Each row contains an individual monitoring observation, with the following schema:

* `country`: Country (Brazil, Indonesia, or Philippines)  
* `site`: 6-letter site code 
* `survey`: Survey (Community Support, Sustainable Fishing Practices, Sustainable Ecosystem, or Sustainable Livelihoods)  
* `before_after`: Binary for before (0) or after (1) the intervention  
* `control_impact`: Binary for control (0) or impact (1) site 
* `replicate`: Replicate for analysis (anonymized indivudal for Community Support, Sustainable Fishing Practices, and Sustainable Livelihoods surveys, and dive site location for Sustainable Ecosystems surveys)  
* `indicator_full`: Full indicator name, breaking community support and sustainable fishing practices indicators out across 6 different behavior changes (Licensing; Catch Reporting; Enforcement; NTZ Compliance; TURF Compliance; Management Participation)  
* `indicator_condensed`: Condensed indicator name, aggregating together community support and sustainable fishing practices indicators across 6 different behavior changes 

### `site_data.csv`  

Each row contains descriptive statistics and demographics for a single site, with the following schema:  

* `country`: Country (Brazil, Indonesia, or Philippines)  
* `site`: 6-letter site code 
* `lat`: Latitude for centroid of site
* `lon`: Longitude for centroid of site
* `site_type`: One of either "Intervention" or "Control"
* `fishers_ff_target_communities`: Number of fishers in Fish Forever target communities  
* `area_current_turfs_ha`: Area of TURFs, as of 2017  
* `area_current_reserves_ha`: Area of reserves (NTZs), as of 2017 


### `survey_question_lookup.csv`

Each row contains the exact survey question that was asked at a particular site, and which survey and indicator the question corresponds to, with the following schema:

* `country`: Country (Brazil, Indonesia, or Philippines)  
* `site_code`: 6-letter site code 
* `survey`: Survey (Community Support, Sustainable Fishing Practices, Sustainable Ecosystem, or Sustainable Livelihoods)  
* `indicator_full`: Full indicator name, breaking community support and sustainable fishing practices indicators out across 6 different behavior changes (Licensing; Catch Reporting; Enforcement; NTZ Compliance; TURF Compliance; Management Participation)  
* `indicator_condensed`: Condensed indicator name, aggregating together community support and sustainable fishing practices indicators across 6 different behavior changes 
* `full_question`: Full survey question, as presented to the survey participant   
* `survey_question_code`: Cleaned version of `full_question` that removes spaces, capital letters, and special characters  
* `response_options`Full survey response options, as presented to the survey participant  

### `phils_matching_data.csv`

Each row contains the attribute score used for base characteristics used during the control site selection for the Philippines sustainable ecosystems and sustainable livelihoods survey, with the following schema:

* `site`: 6-letter site code 
* `control_impact`: Binary for control (0) or impact (1) site 
* `country`: Country (Philippines)  
* `survey`: Survey (Sustainable Ecosystem or Sustainable Livelihoods)  
* `baseline_characteristic`: Name of baseline characteristic  
* `value`: Attribute score (ranked 1 to 5)  

## License

The software code contained within this repository is made available under the [MIT license](http://opensource.org/licenses/mit-license.php). The data and figures are made available under the [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/) license.

**Please note:** To ensure reproducibility and in order to manage package dependencies, we use the `renv` package. When you first clone this repo onto your machine, run `renv::restore()` to ensure you have all correct package versions installed in the project. Please see the `renv` page for more information: https://rstudio.github.io/renv/articles/renv.html.