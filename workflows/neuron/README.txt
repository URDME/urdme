README for worksflows/neuron/

S. Engblom 2019-11-30 (Revision)
S. Engblom 2018-06-28 (Major revision)
A. Senek 2017-05-31

These files form the foundation for simulating a neural membrane as
described in [1] (partially developed in [4,5]). Also distributed are
tools to simulate synaptic transmission from two different model
synapses taken from [2] and described in [4].

This is the first release intended for URDME 1.4. A more mature
interface is scheduled for release with URDME 1.5.

There are minor dependencies on the TREES Toolbox,
www.treestoolbox.org, on STENGLIB, www.stenglib.org, and on Comsol,
www.comsol.com.

Channel modeling, synapse modeling, and spatial modeling are handled
by the functions found within the directories channel/, synapse/, and
spatial/, respectively.

Single modeling examples:
examples/

ion_run
	Investigates an example taken from Rallpack 3.

neuron_run
	Runs a simulation of a realistic neuron using a pre-processed
	tree.

synaptic_run
	Uses a reference voltage from Rallpack 3 to propagate a
	current across a synapse by AMPAa and NMDAa channels, see [3].

Coupled models:
examples/

twoNeurons_run
	Uses a reference solution from Rallpack 3 to propagate a
	potential across a synapse into a neuron, and across a second
	synapse and finally a second neuron.

circleNeurons_run
	An example which can easily be made into a Comsol-model. The
	neuronal structure is made up of 9 straight cable neurons
	connected by 9 synapses at each corner. Using the same
	framework as pyramideTree_run we calculate the currents inside
	the neurons and synapses which we then save as a
	.mat-file. Using Comsol one can then read the currents into
	the simulation, expressed as line current sources, and
	propagate them into space. The solution to the electrical
	field over time can be found in circle_fire.gif, and took
	approx 20 min to calculate.

pyramideTree_run
	An example of a complete coupling of eight neurons. An image
	of the neuronal structure can be found in pyramide_tree.eps,
	and an animation in pyramide_tree.gif. Plots the membrane
	potentials at the end-node of each neuron, see
	pyramide_tree.eps for structure and numbering.

Utility functions:
utils/

model_builder
	Creates necessary parts for discretization of the neuronal
	structure.

neuron2comsol
	Returns a Comsol model in Matlab, to save into an .mph-file
	using, e.g., mphsave(model,'<mph-name>'). The user must ensure
	that Comsol is able to find the function getCurrent on startup
	when running the Comsol solver.

findPredecessor
	Finds the connected vertices for neuron2comsol.

make_gif
	Creates a .gif-animation of a tree with a given dynamic
	potential.
			
ISI_calc
	Calculates statistics (interspike interval), for a firing
	neuron.

References:
  [1] P. Bauer, S. Engblom, S. Mikulovic, and A. Senek: "Multiscale
       modeling via split-step methods in neural firing",
       Math. Comput. Model. Dyn. Syst. 24(4):409--425, (2018).
  [2] A. Destexhe, Z. F. Mainen, and T. J. Sejnowski: "Kinetic models
       of synaptic transmission", Meth. Neur. Model. 2:1--25 (1998).
  [3] A. Destexhe, Z. F. Mainen, and T. J. Sejnowski: "Synthesis of
       models for excitable membranes, synaptic transmission and
       neuromodulation using a common kinetic formalism", J.
       Comput. Neurosci. 1(3):195--230 (1994).
  [4] A. Senek: "Multiscale Stochastic Neuron Modeling: with
       applications in deep brain stimulation", MSc thesis in Engineering
       Physics, Uppsala university (2017).
       http://uu.diva-portal.org/smash/record.jsf?pid=diva2%3A1143600&dswid=-6496
  [5] E. Berwald: "Towards mesoscopic modeling of firing neurons: a
       feasibility study", MSc thesis in Engineering Physics, Uppsala
       university (2014).
       http://uu.diva-portal.org/smash/record.jsf?pid=diva2%3A712486&dswid=-6336
