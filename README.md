# Indoor Temperature Estimation Using Outside Temperature and Class Schedules
We predict the indoor temperature in a building by using the outdoor tempera-
ture through one temperature sensor, which is located outside the building, an
initial indoor temperature, and a normal distribution of students arrival to the
building using their class schedule. We also have an indoor temperature sen-
sor, which is used to collect the ground truth. Our model predicts the indoor
temperature based on when the entrance door to the building is opened and
for how long. To predict the temperature, the model estimates the amount of
air that escaped (based on the temperature difference between the indoor and
outdoor and the time the door was opened). The estimation of the air escaped
is used to compute the indoor temperature after the air escaped. The idea is to
have a normal distribution for each class that predicts the number of students
arriving to class, sum the number of students and coming to the building one
after the other to get the number of seconds the door is opened every minute
and predict the indoor temperature using our model. We use linear regression
to increase the accuracy. We also use linear and quadratic b-splines to predict
the indoor temperature. We compare all these methods to the ground truth.

