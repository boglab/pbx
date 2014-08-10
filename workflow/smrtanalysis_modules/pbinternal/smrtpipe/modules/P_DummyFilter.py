import os
import SMRTpipe.modules.P_Filter as P_Filter_Module
from SMRTpipe.modules.P_Module import P_Module
from SMRTpipe.engine.common import PRE_WORKFLOW
from SMRTpipe.engine.SmrtPipeFiles import cmdLineInput
from SMRTpipe.engine.SmrtPipeTasks import task

# Patches the SMRTpipe P_Filter module to create links to existing files instead of doing filtering again
# This should be placed first in the protocol
# Wish there was a nicer way to do this

@task(inputs={'plsFofn': P_Filter_Module.inputPlsFofn},
      outputs={'rgnFofn': P_Filter_Module.filteredRgnFofn,
               'summary': P_Filter_Module.filterSummary,
               'rgnDir': P_Filter_Module.filteredRgnDir},
      errors={"Unable to locate BaseCalls within hdf5 file":
                  "No Base Calls found in Primary; Pulse2Base may have failed."},
      localOnly=True,
      coreTask=True)
def dummy_filter(self, files):
    
    origin = self.setting("linkOrigin", None)
    
    cmds = []
    
    cmds.append("ln -sf %s %s" % (os.path.join(origin, files.rgnFofn.relativePath), files.rgnFofn.path))
    cmds.append("ln -sf %s %s" % (os.path.join(origin, files.summary.relativePath), files.summary.path))
    cmds.append("ln -sf %s %s" % (os.path.join(origin, files.rgnDir.relativePath), files.rgnDir.path))
    
    return cmds

@task(inputs={'plsFofn': P_Filter_Module.inputPlsFofn, 'rgnFofn': P_Filter_Module.filteredRgnFofn},
      outputs={'subreads': P_Filter_Module.filteredSubreads,
               'subreadFastq': P_Filter_Module.filteredSubreadFastq},
      errors={"Assertion `nDims == 1' failed": "AG4C Primary Analysis Reports are no longer supported."},
      localOnly=True)
def dummy_subreads(self, files):
    
    origin = self.setting("linkOrigin", None)
    
    cmds = []
    
    cmds.append("ln -sf %s %s" % (os.path.join(origin, files.subreads.relativePath), files.subreads.path))
    cmds.append("ln -sf %s %s" % (os.path.join(origin, files.subreadFastq.relativePath), files.subreadFastq.path))
    
    return cmds

@task(inputs={'plsFofn': P_Filter_Module.inputPlsFofn, 'rgnFofn': P_Filter_Module.filteredRgnFofn},
      outputs={'CCSsubreads': P_Filter_Module.filteredCCSSubreads,
               'CCSsubreadFastq': P_Filter_Module.filteredCCSSubreadFastq},
      errors={"Assertion `nDims == 1' failed": "AG4C Primary Analysis Reports are no longer supported."},
      enabled=False,
      toggleable=True,
      localOnly=True
)
def dummy_ccsSubreads(self, files):
    
    origin = self.setting("linkOrigin", None)
    
    cmds = []
    
    cmds.append("ln -sf %s %s" % (os.path.join(origin, files.CCSsubreads.relativePath), files.CCSsubreads.path))
    cmds.append("ln -sf %s %s" % (os.path.join(origin, files.CCSsubreadFastq.relativePath), files.CCSsubreadFastq.path))
    
    return cmds

@task(inputs={'rgnFofn': P_Filter_Module.filteredRgnFofn},
      outputs={'subreadSummary': P_Filter_Module.filterSubreadSummaryCsv},
      localOnly=True)
def dummy_subreadSummary(self, files):
    
    origin = self.setting("linkOrigin", None)
    
    cmds = []
    
    cmds.append("ln -sf %s %s" % (os.path.join(origin, files.subreadSummary.relativePath), files.subreadSummary.path))
    
    return cmds

P_Filter_Module.P_Filter.filter = dummy_filter
P_Filter_Module.P_Filter.subreads = dummy_subreads
P_Filter_Module.P_Filter.ccsSubreads = dummy_ccsSubreads
P_Filter_Module.P_Filter.subreadSummary = dummy_subreadSummary

class P_DummyFilter( P_Module ):

    def validateSettings( self ):
        """Extract relevant settings from the context and store as attributes
        of this module, setting defaults and validating where necessary."""
        errors = P_Module.validateSettings( self )
        return errors

    @task( inputs   = { 'cmdline'   : cmdLineInput },
           when     = PRE_WORKFLOW )
    def runScripts( self, files ):
        return []
