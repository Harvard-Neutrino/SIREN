# Monkey patch DarkNews logger to hide printouts
import functools

dn_has_modelcontainer_logger = False
try:
    from DarkNews.ModelContainer import ModelContainer
    ModelContainer_configure_logger = ModelContainer.configure_logger
    dn_has_modelcontainer_logger = True
except:
    pass

if dn_has_modelcontainer_logger:
    @functools.wraps(ModelContainer.configure_logger)
    def suppress_info(self, logger, loglevel="INFO", prettyprinter=None, logfile=None, verbose=False):
        return ModelContainer_configure_logger(self, logger, loglevel="WARNING", prettyprinter=prettyprinter, logfile=logfile, verbose=verbose)

    ModelContainer.configure_logger = suppress_info

