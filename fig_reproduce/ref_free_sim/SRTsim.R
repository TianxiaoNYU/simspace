library('SRTsim')
save_path <- './'

## Since SRTsim reference-free simulation function is implemented in shiny app.
## and everything has to be done by manual operations in the shiny app,
## we provide the code here to show how to get to the shiny app and save the simulation result.
## Please run the code below to launch the shiny app and follow the instructions in the app
## to finish the simulation and save the result.

## Once within the shiny app, one can always view our generated reference-free simulation
## by selecting the "Load Example Data" option in the "Choose CSV Files" section.
## Through that one can know how we generated the reference-free simulation used in the paper.
shinySRT1 <- SRTsim_shiny()

spatial_res <- shinySRT1[['simInfo']]
write.csv(spatial_res, paste0(save_path, 'ref_free_0.csv'), quote = F)
