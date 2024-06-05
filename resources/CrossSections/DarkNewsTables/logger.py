# Monkey patch DarkNews logger to hide printouts

from DarkNews.ModelContainer import ModelContainer
ModelContainer_configure_logger = ModelContainer.configure_logger

@functools.wraps(ModelContainer.configure_logger)
def suppress_info(self, logger, loglevel="INFO", prettyprinter=None, logfile=None, verbose=False):
    return ModelContainer_configure_logger(self, logger, loglevel="WARNING", prettyprinter=prettyprinter, logfile=logfile, verbose=verbose)

ModelContainer.configure_logger = suppress_info

