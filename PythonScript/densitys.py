import CoolProp.CoolProp as CP

ρ1 = CP.PropsSI("D", "T", 85, "P", (6894.76*370), "O2")
ρ2 = CP.PropsSI("D", "T", 298, "P", (6894.76*370), "Ethanol")
ρ3 = CP.PropsSI("D", "T", 298, "P", (6894.78*2200), "Nitrogen")

print("LOX at 370: "+str(ρ1)+", Ethanol at 370: "+str(ρ2)+", Nitrogen at 2200: "+str(ρ3))
