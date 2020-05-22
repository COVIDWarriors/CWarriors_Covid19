# COVIDWarriors

Current global pandemic due to the SARS-COVID19 have struggled and pushed the limits of global health systems. Supply chain disruptions and scarce availability of RNA extraction kits have motivated worldwide actors to do research on virus detection systems.

The NGO ***COVIDWarriors*** made a donation to multiple hospitals in Spain to automate the PCR testing of SARS-cov2 by suppplying them with opensource programmable robots, developed by the US company ***Opentrons***. They also collaborated in the installation and setup by providing one engineer.

The following methodology prepared in this repository will take advantage of commercial kits to perform the RNA extraction.

The main contributors have developed some multipurpouse coding methodology which has been very useful to prepare the extraction protocols.

--------------
# Coding methodology


- **Classes:** reactive classes have been proposed which include the pipetting characteristics and other parameters related to the labware in which they are placed.

- **Functions:**

  - **Divide volume:** generates a list of pipetting volumes when the volume to be transferred is bigger than the tip capacity.
  - **Custom mix:** generates a multi pipetting function in one same well with the number of repetitions required.
  - **Move volume multichannel:** moves volumes of liquid from one source to a destination with upgraded parameters as a prior rinse to the aspiration, blows out at the end, pipette position (X,Y) dynamic modification and touch tip if requested.
  - **Distribute custom:** distributes a certain volume of reactive within multiple wells with air gap parameters and disposal selection.
  - **Calculate height:** calculates the height from which the pipette must aspirate the reactive taking into account the remaining volume in the source well, by minimizing the tip wetting to avoid droplets. At the same time, if no volume is left in the tube, it will move its sourcing position to the next well defined as a source.

- **General code structure:** coding has been structured as in the protocol diagrams by splitting in very concise steps, fed by the previously defined functions, and controlled by a dictionary type variable which will activate or deactivate the tasks, easing the debugging and fine tunning process of the robot.

--------------
# Robot operation description

One working unit is composed by 4 multidispensing robots divided in 3 different stations named A, B, and C. Each of the stations performs different actions as described below:

- **Station A**: This station does a sample setup. Original samples are distributed in 4 racks of 24 samples each. The samples are redistributed in a 96 deepwell plate and a control reactive is added to each well.

- **Station B**: RNA extraction procedure using magnetic microbeads and transfer to a 96 well elution plate .

- **Station C**: The qPCR plate is prepared by adding the required volume of elution from the elution plate coming from station B and the required volume of Mastermix.

--------------
A truly sincere recognition for their time, support and contribution to:

- Jose Luis Villanueva
- Eva Gonzalez
- Joan Anton Puig

And all the team of the CDB

Not last, thanks to the support to the whole team of COVIDWarriors, with special thanks to:
- Andreu Veà
- Rocío Martínez
- Alex Gasulla
- Ramón Martínez
- Aitor Gastaminza
- Ernesto Sánchez

Disclaimer:

For research use only.
