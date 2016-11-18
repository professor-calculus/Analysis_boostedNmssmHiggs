import os

class SubmitHtcondorExeJob:

	Executable
	Arguments
	Output
	Log
	Error
	Input
	RequestMem	

	def submit():
		submitString = "condor_submit htcondor_executable.job "

		submitString += "-append 'Executable = " + Executable + "' "
		submitString += "-append 'arguments = " + Arguments + "' "
		submitString += "-append 'Output = " + Output + "' "
		submitString += "-append 'Log = " + Log + "' "
		submitString += "-append 'Error = " + Error + "' "
		submitString += "-append 'input = " + Input + "' "
		submitString += "-append 'request_memory = " + RequestMem + "' "

		print submitString
		# os.system("%s" % submitString)



		# want to save this in the dice directory