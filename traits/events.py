from enthought.traits.api import HasTraits, Float, on_trait_change, List, Str

class Person(HasTraits):
    weight_kg = Float
    height_m =  Str
    bmi = Str


    @on_trait_change( 'weight_kg' )
    def grr(self,**sd):
      print "======Weigth was updated"

bob = Person()

bob.configure_traits()
