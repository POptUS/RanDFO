class MissingTrueFunctionError(Exception):
    """Raised when f_true is required for stochastic optimization but not provided."""
    def __init__(self, message="Stochastic mode requires both a noisy function and the true function (f_true)."):
        self.message = message
        super().__init__(self.message)

class ConflictingTrackingOptionsError(Exception):
    """Raised when both track_objfun and track_progress are set to True."""
    def __init__(self, message="Only one of 'track_objfun' or 'track_progress' can be True at a time for better console output"):
        self.message = message
        super().__init__(self.message)

