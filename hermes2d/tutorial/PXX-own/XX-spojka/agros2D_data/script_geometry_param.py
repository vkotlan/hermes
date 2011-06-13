# parametry
r_osa=0.005
r_osa_izolace=0.007
r_magnet_vnitrni=0.015
r_vzduch_1=0.016
r_Fe_rozt=0.025
r_vzduch_2=0.0255
r_Fe_nerozt=0.035
r_vzduch_3=0.037
r_magnet_vnejsi_obal=0.052
r_vzduch_4=0.3
# parametry magnety - geometrie
n=8. 				# magnety tvori n-uhelnik
delta=0.0005		# odsazeni od hran
alpha_n=(n-2)/n*pi 	# vnitrni uhel
alpha_s=2*pi/n		# stredovy uhel
mg_sirka=1./6.*(2.*r_vzduch_3*tan(alpha_n/2))
mg_vyska=(r_magnet_vnejsi_obal*cos(asin(mg_sirka/(2.*r_magnet_vnejsi_obal))))-r_vzduch_3
A_x=r_vzduch_3*cos(alpha_s)
A_y=r_vzduch_3*sin(alpha_s)
dx_B=mg_sirka/2*sin(alpha_s)
dy_B=mg_sirka/2*cos(alpha_s)
dx_C=mg_vyska*cos(alpha_s)
dy_C=mg_vyska*sin(alpha_s)

# startup script
Br = 1.28
mur_mg = 1.11
sigma_Fe_rozt = 6e6
sigma_Fe_nerozt= 2e6
n = 3000
p = 1
f = p*n/60
omega = 2*pi*f

# model
newdocument("spojka", "planar", "magnetic", 1, 2, "disabled", 1, 1, 0, "steadystate", 1, 1, 0)

# boundaries
addboundary("hranice", "magnetic_vector_potential", 0, 0)

# materials
addmaterial("vzduch", 0, 0, 1, 0, 0, 0, 0, 0, 0)
addmaterial("osa", 0, 0, 1, 5e+06, 0, 0, 0, 0, omega)
addmaterial("izolace", 0, 0, 1, 0, 0, 0, 0, 0, omega)
addmaterial("Fe_rozt", 0, 0, 2000, sigma_Fe_rozt, 0, 0, 0, 0, omega)
addmaterial("Fe_nerozt", 0, 0, 350, sigma_Fe_nerozt, 0, 0, 0, 0, 0)
addmaterial("magnet_vnitrni", 0, 0, 1, 0, 0, 0, 0, 0, 0)
addmaterial("magnet_vnejsi_obal", 0, 0, 1, 0, 0, 0, 0, 0, 0)
addmaterial("M_N", 0, 0, mur_mg, 0, Br, 90, 0, 0, 0)
addmaterial("M_NW", 0, 0, mur_mg, 0, Br, 45, 0, 0, 0)
addmaterial("M_W", 0, 0, mur_mg, 0, Br, 0, 0, 0, 0)
addmaterial("M_SW", 0, 0, mur_mg, 0, Br, 315, 0, 0, 0)
addmaterial("M_S", 0, 0, mur_mg, 0, Br, 270, 0, 0, 0)
addmaterial("M_SE", 0, 0, mur_mg, 0, Br, 225, 0, 0, 0)
addmaterial("M_E", 0, 0, mur_mg, 0, Br, 180, 0, 0, 0)
addmaterial("M_NE", 0, 0, mur_mg, 0, Br, 135, 0, 0, 0)

# edges
addedge(0, -r_osa, 0, r_osa, 180, "none")
addedge(0, r_osa, 0, -r_osa, 180, "none")
addedge(0, -r_osa_izolace, 0, r_osa_izolace, 180, "none")
addedge(0, r_osa_izolace, 0, -r_osa_izolace, 180, "none")
addedge(0, -r_magnet_vnitrni, 0, r_magnet_vnitrni, 180, "none")
addedge(0, r_magnet_vnitrni, 0, -r_magnet_vnitrni, 180, "none")
addedge(0, -r_vzduch_1, 0, r_vzduch_1, 180, "none")
addedge(0, r_vzduch_1, 0, -r_vzduch_1, 180, "none")
addedge(0, -r_Fe_rozt, 0, r_Fe_rozt, 180, "none")
addedge(0, r_Fe_rozt, 0, -r_Fe_rozt, 180, "none")
addedge(0, -r_vzduch_2, 0, r_vzduch_2, 180, "none")
addedge(0, r_vzduch_2, 0, -r_vzduch_2, 180, "none")
addedge(0, -r_Fe_nerozt, 0, r_Fe_nerozt, 180, "none")
addedge(0, r_Fe_nerozt, 0, -r_Fe_nerozt, 180, "none")
addedge(0, -r_vzduch_3, 0, r_vzduch_3, 180, "none")
addedge(0, r_vzduch_3, 0, -r_vzduch_3, 180, "none")
addedge(0, -r_magnet_vnejsi_obal, 0, r_magnet_vnejsi_obal, 180, "none")
addedge(0, r_magnet_vnejsi_obal, 0, -r_magnet_vnejsi_obal, 180, "none")
addedge(0, -r_vzduch_4, 0, r_vzduch_4, 180, "hranice")
addedge(0, r_vzduch_4, 0, -r_vzduch_4, 180, "hranice")

addedge(-(mg_sirka/2-delta), r_vzduch_3+delta, (mg_sirka/2-delta), r_vzduch_3+delta, 0, "none")
addedge(-(mg_sirka/2-delta), r_vzduch_3+mg_vyska-delta, (mg_sirka/2-delta), r_vzduch_3+mg_vyska-delta, 0, "none")
addedge(-(mg_sirka/2-delta), r_vzduch_3+delta, -(mg_sirka/2-delta), r_vzduch_3+mg_vyska-delta, 0, "none")
addedge((mg_sirka/2-delta), r_vzduch_3+delta, (mg_sirka/2-delta), r_vzduch_3+mg_vyska-delta, 0, "none")

addedge((A_x-dx_B+sqrt(2)*delta), (A_y+dy_B), (A_x+dx_B), (A_y-dy_B+sqrt(2)*delta), 0, "none")
addedge((A_x-dx_B+dx_C), (A_y+dy_B+dy_C-sqrt(2)*delta), (A_x+dx_B+dx_C-sqrt(2)*delta), (A_y-dy_B+dy_C), 0, "none")
addedge((A_x-dx_B+sqrt(2)*delta), (A_y+dy_B), (A_x-dx_B+dx_C), (A_y+dy_B+dy_C-sqrt(2)*delta), 0, "none")
addedge((A_x+dx_B), (A_y-dy_B+sqrt(2)*delta), (A_x+dx_B+dx_C-sqrt(2)*delta), (A_y-dy_B+dy_C), 0, "none")

addedge(r_vzduch_3+delta, (mg_sirka/2-delta), r_vzduch_3+delta, -(mg_sirka/2-delta), 0, "none")
addedge(r_vzduch_3+mg_vyska-delta, -(mg_sirka/2-delta), r_vzduch_3+mg_vyska-delta, (mg_sirka/2-delta), 0, "none")
addedge(r_vzduch_3+delta, (mg_sirka/2-delta), r_vzduch_3+mg_vyska-delta, (mg_sirka/2-delta), 0, "none")
addedge(r_vzduch_3+delta, -(mg_sirka/2-delta), r_vzduch_3+mg_vyska-delta, -(mg_sirka/2-delta), 0, "none")

addedge((A_x-dx_B+sqrt(2)*delta), -(A_y+dy_B), (A_x+dx_B), -(A_y-dy_B+sqrt(2)*delta), 0, "none")
addedge((A_x-dx_B+dx_C), -(A_y+dy_B+dy_C-sqrt(2)*delta), (A_x+dx_B+dx_C-sqrt(2)*delta), -(A_y-dy_B+dy_C), 0, "none")
addedge((A_x-dx_B+sqrt(2)*delta), -(A_y+dy_B), (A_x-dx_B+dx_C), -(A_y+dy_B+dy_C-sqrt(2)*delta), 0, "none")
addedge((A_x+dx_B), -(A_y-dy_B+sqrt(2)*delta), (A_x+dx_B+dx_C-sqrt(2)*delta), -(A_y-dy_B+dy_C), 0, "none")

addedge(-(mg_sirka/2-delta), -(r_vzduch_3+delta), (mg_sirka/2-delta), -(r_vzduch_3+delta), 0, "none")
addedge(-(mg_sirka/2-delta), -(r_vzduch_3+mg_vyska-delta), (mg_sirka/2-delta), -(r_vzduch_3+mg_vyska-delta), 0, "none")
addedge(-(mg_sirka/2-delta), -(r_vzduch_3+delta), -(mg_sirka/2-delta), -(r_vzduch_3+mg_vyska-delta), 0, "none")
addedge((mg_sirka/2-delta), -(r_vzduch_3+delta), (mg_sirka/2-delta), -(r_vzduch_3+mg_vyska-delta), 0, "none")

addedge(-(A_x-dx_B+sqrt(2)*delta), -(A_y+dy_B), -(A_x+dx_B), -(A_y-dy_B+sqrt(2)*delta), 0, "none")
addedge(-(A_x-dx_B+dx_C), -(A_y+dy_B+dy_C-sqrt(2)*delta), -(A_x+dx_B+dx_C-sqrt(2)*delta), -(A_y-dy_B+dy_C), 0, "none")
addedge(-(A_x-dx_B+sqrt(2)*delta), -(A_y+dy_B), -(A_x-dx_B+dx_C), -(A_y+dy_B+dy_C-sqrt(2)*delta), 0, "none")
addedge(-(A_x+dx_B), -(A_y-dy_B+sqrt(2)*delta), -(A_x+dx_B+dx_C-sqrt(2)*delta), -(A_y-dy_B+dy_C), 0, "none")

addedge(-(r_vzduch_3+delta), (mg_sirka/2-delta), -(r_vzduch_3+delta), -(mg_sirka/2-delta), 0, "none")
addedge(-(r_vzduch_3+mg_vyska-delta), -(mg_sirka/2-delta), -(r_vzduch_3+mg_vyska-delta), (mg_sirka/2-delta), 0, "none")
addedge(-(r_vzduch_3+delta), (mg_sirka/2-delta), -(r_vzduch_3+mg_vyska-delta), (mg_sirka/2-delta), 0, "none")
addedge(-(r_vzduch_3+delta), -(mg_sirka/2-delta), -(r_vzduch_3+mg_vyska-delta), -(mg_sirka/2-delta), 0, "none")

addedge(-(A_x-dx_B+sqrt(2)*delta), (A_y+dy_B), -(A_x+dx_B), (A_y-dy_B+sqrt(2)*delta), 0, "none")
addedge(-(A_x-dx_B+dx_C), (A_y+dy_B+dy_C-sqrt(2)*delta), -(A_x+dx_B+dx_C-sqrt(2)*delta), (A_y-dy_B+dy_C), 0, "none")
addedge(-(A_x-dx_B+sqrt(2)*delta), (A_y+dy_B), -(A_x-dx_B+dx_C), (A_y+dy_B+dy_C-sqrt(2)*delta), 0, "none")
addedge(-(A_x+dx_B), (A_y-dy_B+sqrt(2)*delta), -(A_x+dx_B+dx_C-sqrt(2)*delta), (A_y-dy_B+dy_C), 0, "none")

# labels
addlabel(((r_vzduch_4-r_magnet_vnejsi_obal)/2+r_magnet_vnejsi_obal), 0, 0, 0, "vzduch")
addlabel((r_magnet_vnejsi_obal-3*delta), 0, 0, 0, "magnet_vnejsi_obal")
addlabel(((r_vzduch_3-r_Fe_nerozt)/2+r_Fe_nerozt), 0, 0, 0, "vzduch")
addlabel(((r_Fe_nerozt-r_vzduch_2)/2+r_vzduch_2), 0, 0, 0, "Fe_nerozt")
addlabel(((r_vzduch_2-r_Fe_rozt)/2+r_Fe_rozt), 0, 0, 0, "vzduch")
addlabel(((r_Fe_rozt-r_vzduch_1)/2+r_vzduch_1), 0, 0, 0, "Fe_rozt")
addlabel(((r_vzduch_1-r_magnet_vnitrni)/2+r_magnet_vnitrni), 0, 0, 0, "vzduch")
addlabel(((r_magnet_vnitrni-r_osa_izolace)/2+r_osa_izolace), 0, 0, 0, "magnet_vnitrni")
addlabel(((r_osa_izolace-r_osa)/2+r_osa), 0, 0, 0, "izolace")
addlabel(r_osa/2, 0, 0, 0, "osa")
addlabel(0, r_vzduch_3+mg_vyska/2, 0, 0, "M_N")
addlabel((A_x+dx_B/2), (A_y+dy_B/2), 0, 0, "M_NW")
addlabel(r_vzduch_3+mg_vyska/2, 0, 0, 0, "M_W")
addlabel((A_x+dx_B/2), -(A_y+dy_B/2), 0, 0, "M_SW")
addlabel(0, -(r_vzduch_3+mg_vyska/2), 0, 0, "M_S")
addlabel(-(A_x+dx_B/2), -(A_y+dy_B/2), 0, 0, "M_SE")
addlabel(-(r_vzduch_3+mg_vyska/2), 0, 0, 0, "M_E")
addlabel(-(A_x+dx_B/2), (A_y+dy_B/2), 0, 0, "M_NE")