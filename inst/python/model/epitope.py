class BEpitope:
    def __init__(self, start, end, max_fitness):
        self.start = start
        self.end = end
        self.max_fitness = max_fitness

    def __repr__(self):
        return "(%s to %s, maxfit: %s)" % (
            self.start, self.end, self.max_fitness)
            
class TEpitope:
    def __init__(self, start, escape_positions, recognized_aa):
        self.start = start
        self.end = start + max(escape_positions)*3
        self.escape_positions = dict(zip(escape_positions, recognized_aa))
        
    def __repr__(self):
        return "(start: %s, escape positions: %s)" % (self.start, self.escape_positions)
