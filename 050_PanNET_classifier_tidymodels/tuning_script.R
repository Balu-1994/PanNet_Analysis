library(tidymodels)
library(future)

plan(multicore, workers = 8)

NDD_TD = readRDS("/home/rstudio/NDD_TD.Rds")

set.seed = (1618)
NDD_TDfolds = nested_cv(NDD_TD,
                        outside = vfold_cv(v = 4, 
                                           repeats = 5,
                                           strata = CC_Epi_newLRO),
                        inside = bootstraps(time = 10,
                                            strata = CC_Epi_newLRO))

rf_ranger = rand_forest(mode = "classification",
                        trees = tune(),
                        mtry = tune()) %>% 
  set_engine("ranger")

rf_grid = grid_regular(trees(c(100, 500)),
                       mtry(c(25, 125)), 
                       levels = 5)

NDD_recipes = readRDS("/home/rstudio/NDD_recipes.Rds")

DMP_all_ranger = workflow() %>% 
  add_model(rf_ranger) %>% 
  add_recipe(NDD_recipes$DMP_all)

tune_rf = function(object, wf, tg){
  wf %>% 
    tune_grid(object, 
              grid = tg,
              control = control_grid(parallel_over = "everything"))
  
}

#tic()
set.seed(1648)
rf_ranger_tune = map(test_folds$inner_resamples, 
                     tune_rf,
                     wf = test_workflow,
                     tg = test_grid)
#toc()