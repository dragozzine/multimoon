class Test:
    instance = None
    mass = 10
    dog = 7
    
    def getInstance(self):
        if self.instance == None:
            self.instance = Test()
        return self.instance    