import os

class htcondorExeJob:

	def __init__(self):
		self.Executable = []
		self.Arguements = []
		self.Output = []
		self.Log = []
		self.Error = []
		self.RequestMem = []
		self.Input = []
		self.OtherInputFiles = []


	def add_Executable(self, exe):
		self.Executable = exe

	def add_Arguments(self, arg):
		self.Arguments = arg

	def add_Output(self, out):
		self.Output = out

	def add_Log(self, log):
		self.Log = log

	def add_Error(self, err):
		self.Error = err

	def add_Input(self, inp):
		self.Input = inp

	def add_OtherInputFiles(self, otherInp):
		self.OtherInputFiles = otherInp

	def add_RequestMem(self, reqm):
		self.RequestMem = reqm


	def submitJob(self):

		submitString = "condor_submit $CMSSW_BASE/src/Analysis/Analysis_boostedNmssmHiggs/batch/htcondor_executable.job "

		submitString += "-append 'Executable = " + self.Executable + "' "
		submitString += "-append 'arguments = " + self.Arguments + "' "
		submitString += "-append 'Output = " + self.Output + "' "
		submitString += "-append 'Log = " + self.Log + "' "
		submitString += "-append 'Error = " + self.Error + "' "
		submitString += "-append 'request_memory = " + self.RequestMem + "' "
		submitString += "-append 'input = " + self.Input + "' "
		submitString += "-append 'transfer_input_files = " + self.OtherInputFiles + "' "

		# print submitString
		# os.system("%s" % submitString)
		return os.popen("%s" % submitString).read().rstrip()