h_start = 0.287 #[m]

h_slutt = 0.155 #[m]

m_kule = 0.037 #[kg]

V_slutt = [1.391,1.377,1.35828 ,1.36367, 1.407, 1.361, 1.343, 1.345, 1.3435, 1.343, 1.339] #[m/s]

rulletid = [1.55, 1.54, 1.500, 1.5333, 1.500, 1.533, 1.500, 1.567, 1.533, 1.433] #[s]

g = 9.81 #[m/s**2]
def tap_mekanisk_energi(): 
    Energi_sum = 0
    #finner gjennomsnittstap på de ti forsøkene
    
    for i in V_slutt:
        Energi_sum = h_start*g*m_kule - (i**2 * 0.5 * m_kule + m_kule*g*h_slutt)
        
    
    Energitap_gjennomsnitt = Energi_sum/10
    
    print(Energitap_gjennomsnitt)
    
    return Energitap_gjennomsnitt 