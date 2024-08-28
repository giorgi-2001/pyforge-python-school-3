import logging

logger = logging.getLogger("Event logger")


logger.setLevel(logging.INFO)


formatter = logging.Formatter(
    "%(asctime)s\t%(levelname)s\t%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)


file_handler = logging.FileHandler("/app/src/events.log")

file_handler.setFormatter(formatter)


logger.addHandler(file_handler)


class LogFactory:

    def __init__(
        self,
        status_code: int,
        ip_adress: str,
        method: str,
        path: str,
        process_time: float,
        message: str | None = None
    ):
        self._status_code = status_code
        self._ip_adress = ip_adress
        self._method = method
        self._path = path
        self._process_time = process_time
        self._message = message
    
    @property
    def base_log(self):
        return(
            f"{self._status_code}\t{self._ip_adress}\t"
            f"{self._method}\t{self._path}\t{round(self._process_time, 3)}"
        )
    
    @property
    def message_log(self):
        return (
            f"{self._status_code}\t{self._ip_adress}\t"
            f"{self._method}\t{self._path}\t{round(self._process_time, 3)}\n"
            f"\t/-- {self._message} --/\n"
        )
    
    def set_message(self, message: str):
        self._message = message
    
