source("tree_simulation/model_setup.r")
source("tree_simulation/model_simulate.r")


model1 <- coale_model(20, 1000)
model1 <- model_simulator(model1)

model2 <- coale_model(20, 1000)
model2 <- add_recombination(model2, 0.001)
model2 <- model_simulator(model2)

model3 <- coale_model(20, 1000)
model3 <- add_mutation(model3, 0.001)
model3 <- add_recombination(model3, 0.001)
model3 <- model_simulator(model3)

model4 <- add_bacteria_recomb(model1, 0.001, 0.2)
model4 <- bacteria_recomb_simulator(model4)
