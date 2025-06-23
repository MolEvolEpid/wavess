from model.virus import HIV

class InfectedCD4:
    def __init__(self, infecting_virus: HIV, active: bool):
        assert isinstance(
            active, bool
        ), "Provide bool value for whether CD4 T cell is active (vs latent)"
        self.active = active
        self.infecting_virus = infecting_virus

    def __repr__(self):
        return "Infected CD4. Active: %s. %s" % (
            str(self.active),
            str(self.infecting_virus),
        )

    def become_latent(self):
        self.active = False  # Set to False to indicate the cell is latent

    def become_active(self):
        self.active = True  # Set to True to indicate the cell is actively proliferating
