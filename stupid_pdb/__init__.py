import logging
import sys

# Set up a basic console logger.
# All modules in the package can then get this logger instance.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] [%(name)s] %(message)s",
    stream=sys.stdout,
)

# Get the root logger for this package
logger = logging.getLogger(__name__)

logger.debug("stupid_pdb package initialized.")
