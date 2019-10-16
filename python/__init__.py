from icecube.load_pybindings import load_pybindings
import icecube.icetray, icecube.dataclasses, icecube.phys_services, icecube.tableio
load_pybindings(__name__,__path__)

"""class SSplitter(icecube.icetray.I3Module):
	def S(self,frame):
		self.PushFrame(frame)
		dummyFrame=icecube.icetray.I3Frame(icecube.icetray.I3Frame.Physics)
		header=icecube.dataclasses.I3EventHeader()
		header.event_id=self.count
		self.count+=1
		header.sub_event_stream="SimulationSettings"
		dummyFrame.Put("I3EventHeader",header)
		self.PushFrame(dummyFrame)
	
	def __init__(self,context):
		icecube.icetray.I3Module.__init__(self,context)
		self.Register(icecube.icetray.I3Frame.Stream('S'),self.S)
		self.count=0
		self.AddOutBox("OutBox")
	
	def Configure(self):
		pass

@icecube.icetray.traysegment
def WriteInjectorProperties(tray,name,Filename="InjectorProperties.h5"):
	tray.AddModule(SSplitter)
	from icecube import tableio, hdfwriter, LeptonInjector
	tray.AddModule(tableio.I3TableWriter, name+"_hdfwriter")(
		("SubEventStreams",["SimulationSettings"]),
		("tableservice",hdfwriter.I3HDFTableService(Filename)),
		("keys",["LeptonInjectorProperties"])
	)
	def CleanUp(frame):
		if(not frame.Has("I3EventHeader")):
			return True
		return frame["I3EventHeader"].sub_event_stream!="SimulationSettings"
	tray.AddModule(CleanUp)
"""

"""
#Since LeptonInjector creates only one type of final state at a time, these segments 
#are intended to provide an easy way to sample all of the final states for the 
#interactions of a given neutrino flavor.
#Note that the fractions of events which go into each final state need not be exactly
#physical, since the weighting to a physical spectrum will fill in the details

@icetray.traysegment
def CompleteNuMu(tray,totalEvents,**kwargs):
	ccEvents=int(3.*totalEvents/4.)
	ncEvents=totalEvents-ccEvents
	tray.AddModule(RangedLeptonInjector,
	               NEvents=ccEvents,
	               FinalType1=dataclasses.I3Particle.MuMinus,
	               FinalType2=dataclasses.I3Particle.Hadrons,
	               UseAntiParticles=True,
	               #CrossSectionFile=???,
	               **kwargs)
	
	tray.AddModule(VolumeLeptonInjector,
	               NEvents=ncEvents,
	               FinalType1=dataclasses.I3Particle.Hadrons,
	               FinalType2=dataclasses.I3Particle.NuMu,
	               UseAntiParticles=True,
	               #CrossSectionFile=???,
	               **kwargs)

@icetray.traysegment
def CompleteNuE(tray,totalEvents,includeGlashow=True,**kwargs):
	gr_states = [
			     (dataclasses.I3Particle.Hadrons,dataclasses.I3Particle.Hadrons),
			     (dataclasses.I3Particle.EMinus,dataclasses.I3Particle.NuEBar),
			     (dataclasses.I3Particle.MuMinus,dataclasses.I3Particle.NuMuBar),
			     (dataclasses.I3Particle.TauMinus,dataclasses.I3Particle.NuTauBar)
			    ]
	
	grEvents=0
	if includeGlashow:
		grEvents=int(totalEvents/10.)
	ccEvents=int(3.*(totalEvents-grEvents)/4.)
	ncEvents=totalEvents-grEvents-ccEvents
	tray.AddModule(VolumeLeptonInjector,
	               NEvents=ccEvents,
	               FinalType1=dataclasses.I3Particle.EMinus,
	               FinalType2=dataclasses.I3Particle.Hadrons,
	               UseAntiParticles=True,
	               #CrossSectionFile=???,
	               **kwargs)
	
	tray.AddModule(VolumeLeptonInjector,
	               NEvents=ncEvents,
	               FinalType1=dataclasses.I3Particle.Hadrons,
	               FinalType2=dataclasses.I3Particle.NuE,
	               UseAntiParticles=True,
	               #CrossSectionFile=???,
	               **kwargs)
	
	if includeGlashow:
		for state in gr_states:
			tray.AddModule(VolumeLeptonInjector,
			               NEvents=grEvents/len(gr_states),
			               FinalType1=state[0],
			               FinalType2=state[1],
			               UseAntiParticles=True,
			               #CrossSectionFile=???,
			               **kwargs)
	
@icetray.traysegment
def CompleteNuTau(tray,totalEvents,**kwargs):
	ccEvents=int(3.*totalEvents/4.)
	ncEvents=totalEvents-ccEvents
	tray.AddModule(RangedLeptonInjector,
	               NEvents=ccEvents,
	               FinalType1=dataclasses.I3Particle.TauMinus,
	               FinalType2=dataclasses.I3Particle.Hadrons,
	               UseAntiParticles=True,
	               #CrossSectionFile=???,
	               **kwargs)
	
	tray. (VolumeLeptonInjector,
	               NEvents=ncEvents,
	               FinalType1=dataclasses.I3Particle.Hadrons,
	               FinalType2=dataclasses.I3Particle.NuTau,
	               UseAntiParticles=True,
	               #CrossSectionFile=???,
	               **kwargs)

@icetray.traysegment
def CompleteNeutrino(tray,totalEvents,**kwargs):
	nuEEvents=int(totalEvents/3)
	nuMuEvents=nuEEvents
	nuTauEvents=totalEvents-nuEEvents-nuMuEvents
	CompleteNuE(tray,nuEEvents,includeGlashow=True,**kwargs)
	CompleteNuMu(tray,nuEEvents,**kwargs)
	CompleteNuTau(tray,nuEEvents,**kwargs)
"""
