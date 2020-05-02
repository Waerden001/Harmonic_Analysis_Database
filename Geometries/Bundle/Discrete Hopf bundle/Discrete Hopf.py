class bundle():
    """
    A class to handle discrete bundles with connections and gauge symmetry
    """
    def __init__(self, total_space,base_space, morphism):
        self.total_space = total_space
        self.base_space = base_space
        self.morphism = morphism
        
    def curvature(self):
        pass

    def connection(self):
        pass

    def chern(self):
        pass