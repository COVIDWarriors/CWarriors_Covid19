import math
from opentrons.types import Point
from opentrons import protocol_api
from opentrons.drivers.rpi_drivers import gpio
import time
import os
import numpy as np
from timeit import default_timer as timer
import json
from datetime import datetime
import csv

# metadata
metadata = {
    'protocolName': 'S2 Station B Version 4',
    'author': 'Aitor Gastaminza & José Luis Villanueva & Alex Gasulla',
    'source': 'Hospital Clínic Barcelona & HU Vall Hebrón',
    'apiLevel': '2.1',
    'description': 'Protocol for RNA extraction'
}

NUM_SAMPLES = 8
sample_volume = 200 # Sample volume received in station A
set_temp_on = False # Do you want to start temperature module?
recycle_tip = False # Do you want to recycle tips? It shoud only be set True for testing

#mag_height = 11 # Height needed for NUNC deepwell in magnetic deck
mag_height = 14 # Height needed for NEST deepwell in magnetic deck
temperature = 23


L_deepwell = 8 # Deepwell lenght (NEST deepwell)
#D_deepwell = 8.35 # Deepwell diameter (NUNC deepwell)
multi_well_rack_area = 8 * 71 #Cross section of the 12 well reservoir
deepwell_cross_section_area = L_deepwell ** 2 # deepwell square cross secion area

num_cols = math.ceil(NUM_SAMPLES / 8) # Columns we are working on

def run(ctx: protocol_api.ProtocolContext):

    #Change light to red
    gpio.set_button_light(1,0,0)

    ctx.comment('Actual used columns: '+str(num_cols))
    STEP = 0
    STEPS = { #Dictionary with STEP activation, description, and times
            1:{'Execute': False, 'description': 'Mix beads'},# REMOVE
            2:{'Execute': True, 'description': 'Transfer lysis'},#
            3:{'Execute': True, 'description': 'Wait with magnet OFF', 'wait_time': 900}, #300
            4:{'Execute': True, 'description': 'Incubate wait with magnet ON', 'wait_time': 300}, #300
            5:{'Execute': True, 'description': 'Remove supernatant'},#
            6:{'Execute': True, 'description': 'Switch off magnet'},#
            7:{'Execute': True, 'description': 'Add VHB/WB1'},#
            8:{'Execute': True, 'description': 'Incubate wait with magnet ON', 'wait_time': 300},#300
            9:{'Execute': True, 'description': 'Remove supernatant'},#
            10:{'Execute': True, 'description': 'Switch off magnet'},#
            11:{'Execute': True, 'description': 'Add SPR/WB2'},#
            12:{'Execute': True, 'description': 'Incubate wait with magnet ON', 'wait_time': 300},#300
            13:{'Execute': True, 'description': 'Remove supernatant'},#
            14:{'Execute': True, 'description': 'Switch off magnet'},#
            15:{'Execute': True, 'description': 'Add SPR/WB2'},#
            16:{'Execute': True, 'description': 'Incubate wait with magnet ON', 'wait_time': 300},#300
            17:{'Execute': True, 'description': 'Remove supernatant'},#
            18:{'Execute': True, 'description': 'Allow to dry', 'wait_time': 900},#900
            19:{'Execute': True, 'description': 'Switch off magnet'},#
            20:{'Execute': True, 'description': 'Add water'},#
            21:{'Execute': True, 'description': 'Wait with magnet OFF', 'wait_time': 300},#600
            22:{'Execute': True, 'description': 'Incubate wait with magnet ON', 'wait_time': 300},#300
            23:{'Execute': True, 'description': 'Transfer to final elution plate'},
            }

    """if not ctx.is_simulating():
        folder_path='/data/log_times/'
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path=folder_path+'/time_log.json'"""

    #Define Reagents as objects with their properties
    class Reagent:
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense, flow_rate_aspirate_mix, flow_rate_dispense_mix,
        air_gap_vol_bottom, air_gap_vol_top, disposal_volume, rinse, max_volume_allowed, reagent_volume, reagent_reservoir_volume, num_wells, h_cono, v_fondo, tip_recycling = 'none'):
            self.name = name
            self.flow_rate_aspirate = flow_rate_aspirate
            self.flow_rate_dispense = flow_rate_dispense
            self.flow_rate_aspirate_mix = flow_rate_aspirate_mix
            self.flow_rate_dispense_mix = flow_rate_dispense_mix
            self.air_gap_vol_bottom = air_gap_vol_bottom
            self.air_gap_vol_top = air_gap_vol_top
            self.disposal_volume = disposal_volume
            self.rinse = bool(rinse)
            self.max_volume_allowed = max_volume_allowed
            self.reagent_volume = reagent_volume
            self.reagent_reservoir_volume = reagent_reservoir_volume
            self.num_wells = num_wells
            self.col = 0
            self.vol_well = 0
            self.h_cono = h_cono
            self.v_cono = v_fondo
            self.tip_recycling = tip_recycling
            self.vol_well_original = reagent_reservoir_volume / num_wells

    #Reagents and their characteristics
    Lysis = Reagent(name = 'Lysis',
                    flow_rate_aspirate = 3, # Original = 0.5
                    flow_rate_dispense = 3, # Original = 1
                    flow_rate_aspirate_mix = 15, # Original = 1.5
                    flow_rate_dispense_mix = 25,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = True,
                    max_volume_allowed = 180,
                    reagent_volume = 530, # reagent volume needed per sample
                    reagent_reservoir_volume =  (NUM_SAMPLES + 5) * 530, #70000, #51648
                    num_wells = math.ceil((NUM_SAMPLES + 5) * 530 / 13000), #num_Wells max is 4, 13000 is the reservoir max volume (eventhough reservoir allows 15000)
                    h_cono = 1.95,
                    v_fondo = 750, #1.95 * multi_well_rack_area / 2, #Prismatic
                    tip_recycling = 'A1')

    VHB = Reagent(name = 'VHB',
                    flow_rate_aspirate = 3,
                    flow_rate_dispense = 3,
                    flow_rate_aspirate_mix = 15,
                    flow_rate_dispense_mix = 25,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = True,
                    max_volume_allowed = 180,
                    reagent_volume = 350,
                    reagent_reservoir_volume = (NUM_SAMPLES + 5) * 350, #60000, #38400
                    num_wells = math.ceil((NUM_SAMPLES + 5) * 350 / 13000), #num_Wells max is 4
                    h_cono = 1.95,
                    v_fondo = 750, #1.95 * multi_well_rack_area / 2, #Prismatic
                    tip_recycling = 'A1')

    Beads_PK = Reagent(name = 'Magnetic beads+PK',
                    flow_rate_aspirate = 1,
                    flow_rate_dispense = 1.5,
                    flow_rate_aspirate_mix = 1.5,
                    flow_rate_dispense_mix = 5,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = True,
                    max_volume_allowed = 180,
                    reagent_volume = 500,
                    reagent_reservoir_volume = NUM_SAMPLES * 500, #11920,
                    num_wells = math.ceil((NUM_SAMPLES + 5) * 500 / 13000), #num_Wells max is 4,
                    h_cono = 1.95,
                    v_fondo = 750, #1.95 * multi_well_rack_area / 2, #Prismatic
                    tip_recycling = 'A2')

    SPR = Reagent(name = 'SPR',
                    flow_rate_aspirate = 3, # Original = 1
                    flow_rate_dispense = 3, # Original = 1
                    flow_rate_aspirate_mix = 15,
                    flow_rate_dispense_mix = 25,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = True,
                    max_volume_allowed = 180,
                    reagent_volume = 350,
                    reagent_reservoir_volume = (NUM_SAMPLES + 5) * 350, #120000, #96000
                    num_wells = math.ceil((NUM_SAMPLES + 5) * 350 / 13000), #num_Wells max is 4
                    h_cono = 1.95,
                    v_fondo = 750, #1.95 * multi_well_rack_area / 2, #Prismatic
                    tip_recycling = 'A3')

    Water = Reagent(name = 'Water',
                    flow_rate_aspirate = 3,
                    flow_rate_dispense = 3,
                    flow_rate_aspirate_mix = 15,
                    flow_rate_dispense_mix = 25,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = False,
                    max_volume_allowed = 150,
                    reagent_volume = 50,
                    reagent_reservoir_volume = (NUM_SAMPLES + 5) * 50,
                    num_wells = 1, #math.ceil((NUM_SAMPLES + 5) * 50 / 13000), #num_Wells max is 1
                    h_cono = 1.95,
                    v_fondo = 750) #1.95*multi_well_rack_area/2) #Prismatic

    Elution = Reagent(name = 'Elution',
                    flow_rate_aspirate = 3, # Original 0.5
                    flow_rate_dispense = 3, # Original 1
                    flow_rate_aspirate_mix = 15,
                    flow_rate_dispense_mix = 25,
                    air_gap_vol_bottom = 5,
                    air_gap_vol_top = 0,
                    disposal_volume = 1,
                    rinse = False,
                    max_volume_allowed = 150,
                    reagent_volume = 50,
                    reagent_reservoir_volume = (NUM_SAMPLES + 5) * 50, #14800,
                    num_wells = num_cols, #num_cols comes from available columns
                    h_cono = 4,
                    v_fondo = 4 * math.pi * 4**3 / 3) #Sphere

    Lysis.vol_well = Lysis.vol_well_original
    VHB.vol_well = VHB.vol_well_original
    #Beads_PK.vol_well = Beads_PK.vol_well_original
    SPR.vol_well = SPR.vol_well_original
    Water.vol_well = Water.vol_well_original
    Elution.vol_well = 350 # Arbitrary value

    ctx.comment(' ')
    ctx.comment('###############################################')
    ctx.comment('VOLUMES FOR ' + str(NUM_SAMPLES) + ' SAMPLES')
    ctx.comment(' ')
    ctx.comment('Lysis: ' + str(Lysis.num_wells) + ' wells from well 1 in reservoir 1 with volume ' + str(Lysis.vol_well_original) + ' uL each one')
    ctx.comment('WH1: ' + str(VHB.num_wells) + ' wells from well 5 in reservoir 1 with volume ' + str(VHB.vol_well_original) + ' uL each one')
    ctx.comment('WH2: ' + str(SPR.num_wells) + ' wells from well 9 in reservoir 1 with volume ' + str(SPR.vol_well_original) + ' uL each one')
    ctx.comment('Water: ' + str(Water.num_wells) + ' wells from well 1 in reservoir 2 with volume ' + str(Water.vol_well_original) + ' uL each one')
    ctx.comment('###############################################')
    ctx.comment(' ')

    ###################
    #Custom functions
    def custom_mix(pipet, reagent, location, vol, rounds, blow_out, mix_height, offset):
        '''
        Function for mix in the same location a certain number of rounds. Blow out optional. Offset
        can set to 0 or a higher/lower value which indicates the lateral movement
        '''
        if mix_height == 0:
            mix_height = 1
        pipet.aspirate(1, location = location.bottom(z = mix_height), rate = reagent.flow_rate_aspirate_mix)
        for _ in range(rounds):
            pipet.aspirate(vol, location = location.bottom(z = mix_height), rate = reagent.flow_rate_aspirate_mix)
            pipet.dispense(vol, location = location.top(z = -5).move(Point(x = offset)), rate = reagent.flow_rate_dispense_mix)
        pipet.dispense(1, location = location.bottom(z = mix_height), rate = reagent.flow_rate_dispense_mix)
        if blow_out == True:
            pipet.blow_out(location.top(z = -2)) # Blow out

    def calc_height(reagent, cross_section_area, aspirate_volume):
        nonlocal ctx
        ctx.comment('Remaining volume ' + str(reagent.vol_well) +
                    '< needed volume ' + str(aspirate_volume) + '?')
        if reagent.vol_well < aspirate_volume:
            ctx.comment('Next column should be picked')
            ctx.comment('Previous to change: ' + str(reagent.col))
            # column selector position; intialize to required number
            reagent.col = reagent.col + 1
            ctx.comment(str('After change: ' + str(reagent.col)))
            reagent.vol_well = reagent.vol_well_original
            ctx.comment('New volume:' + str(reagent.vol_well))
            height = (reagent.vol_well - aspirate_volume - reagent.v_cono) / cross_section_area
                    #- reagent.h_cono
            reagent.vol_well = reagent.vol_well - aspirate_volume
            ctx.comment('Remaining volume:' + str(reagent.vol_well))
            if height < 5:
                height = 1
            col_change = True
        else:
            height = (reagent.vol_well - aspirate_volume - reagent.v_cono) / cross_section_area
            reagent.vol_well = reagent.vol_well - aspirate_volume
            ctx.comment('Calculated height is ' + str(height))
            if height < 5:
                height = 1
            ctx.comment('Used height is ' + str(height))
            col_change = False
        return height, col_change

    def move_vol_multi(pipet, reagent, source, dest, vol, x_offset_source, x_offset_dest, pickup_height, rinse, avoid_droplet, wait_time, blow_out):
        # Rinse before aspirating
        if rinse == True:
            #pipet.aspirate(air_gap_vol_top, location = source.top(z = -5), rate = reagent.flow_rate_aspirate) #air gap
            custom_mix(pipet, reagent, location = source, vol = vol, rounds = 20, blow_out = False, mix_height = 0, offset = 0)
            #pipet.dispense(air_gap_vol_top, location = source.top(z = -5), rate = reagent.flow_rate_dispense)

        # SOURCE
        if reagent.air_gap_vol_top != 0: #If there is air_gap_vol, switch pipette to slow speed
            pipet.move_to(source.top(z = 0))
            pipet.air_gap(reagent.air_gap_vol_top) #air gap
            #pipet.aspirate(reagent.air_gap_vol_top, source.top(z = -5), rate = reagent.flow_rate_aspirate) #air gap

        s = source.bottom(pickup_height).move(Point(x = x_offset_source))
        pipet.aspirate(vol, s) # aspirate liquid

        if reagent.air_gap_vol_bottom != 0: #If there is air_gap_vol, switch pipette to slow speed
            pipet.move_to(source.top(z = 0))
            pipet.air_gap(reagent.air_gap_vol_bottom) #air gap
            #pipet.aspirate(air_gap_vol_bottom, source.top(z = -5), rate = reagent.flow_rate_aspirate) #air gap

        if wait_time != 0:
            ctx.delay(seconds=wait_time, msg='Waiting for ' + str(wait_time) + ' seconds.')

        if avoid_droplet == True: # Touch the liquid surface to avoid droplets
            ctx.comment("Moving to: " + str(pickup_height))
            pipet.move_to(source.bottom(pickup_height))

        # GO TO DESTINATION
        d = dest.top(z = -5).move(Point(x = x_offset_dest))
        pipet.dispense(vol - reagent.disposal_volume + reagent.air_gap_vol_bottom, d, rate = reagent.flow_rate_dispense)

        if wait_time != 0:
            ctx.delay(seconds=wait_time, msg='Waiting for ' + str(wait_time) + ' seconds.')

        if reagent.air_gap_vol_top != 0:
            pipet.dispense(reagent.air_gap_vol_top, dest.top(z = 0), rate = reagent.flow_rate_dispense)

        if blow_out == True:
            pipet.blow_out(dest.top(z = 0))

        #if reagent.air_gap_vol_bottom != 0:
        #    pipet.move_to(dest.top(z = 0))
        #    pipet.air_gap(reagent.air_gap_vol_bottom) #air gap
            #pipet.aspirate(air_gap_vol_bottom, dest.top(z = 0),rate = reagent.flow_rate_aspirate) #air gap

    ##########
    # pick up tip and if there is none left, prompt user for a new rack
    def pick_up(pip):
        nonlocal tip_track
        #if not ctx.is_simulating():
        if tip_track['counts'][pip] >= tip_track['maxes'][pip]:
            ctx.pause('Replace ' + str(pip.max_volume) + 'µl tipracks before \
            resuming.')
            pip.reset_tipracks()
            tip_track['counts'][pip] = 0
        pip.pick_up_tip()

    ##########
    def find_side(col):
        if col%2 == 0:
            side = -1 # left
        else:
            side = 1 # right
        return side

####################################
    # load labware and modules
    ######## 12 well rack
    reagent_res = ctx.load_labware('nest_12_reservoir_15ml', '2','reagent deepwell plate 1')
    reagent_res_2 = ctx.load_labware('nest_12_reservoir_15ml', '3','reagent deepwell plate 2')

############################################
    ########## tempdeck
    tempdeck = ctx.load_module('tempdeck', '1')
    if set_temp_on == True:
        tempdeck.set_temperature(temperature)

##################################
    ####### Elution plate - final plate, goes to C
    elution_plate = tempdeck.load_labware(
        'biorad_96_alum',
        'cooled elution plate')

############################################
    ######## Elution plate - comes from A
    magdeck = ctx.load_module('magdeck', '4')
    #deepwell_plate = magdeck.load_labware('nest_96_wellplate_2000ul', 'NEST 96 Deep Well Plate 2 mL') # Change to NEST deepwell plate
    deepwell_plate = magdeck.load_labware('nest_96_wellplate_2000ul', 'NEST 96 Well Plate 2000 µL') # Change to NEST deepwell plate.
    magdeck.disengage()

####################################
    ######## Waste reservoir
    waste_reservoir = ctx.load_labware('nest_1_reservoir_195ml', '5', 'waste reservoir') # Change to our waste reservoir
    waste = waste_reservoir.wells()[0] # referenced as reservoir

####################################
    ######### Load tip_racks
    tips300 = [ctx.load_labware('opentrons_96_tiprack_300ul', slot, '200µl filter tiprack')
        for slot in ['6', '7', '8', '9', '10', '11']]
    #tips1000 = [ctx.load_labware('opentrons_96_filtertiprack_1000ul', slot, '1000µl filter tiprack')
    #    for slot in ['10']]

###############################################################################
    #Declare which reagents are in each reservoir as well as deepwell and elution plate
    Lysis.reagent_reservoir = reagent_res.rows()[0][:4] # 4 columns
    #Beads_PK.reagent_reservoir = reagent_res.rows()[0][Lysis.num_wells:(Lysis.num_wells+Beads_PK.num_wells)]
    VHB.reagent_reservoir = reagent_res.rows()[0][4:8]
    #SPR.reagent_reservoir = reagent_res.rows()[0][VHB.num_wells:(Lysis.num_wells + Beads_PK.num_wells + VHB.num_wells + SPR.num_wells)]
    SPR.reagent_reservoir = reagent_res_2.rows()[0][:8]
    Water.reagent_reservoir = reagent_res.rows()[0][-1]
    work_destinations = deepwell_plate.rows()[0][:Elution.num_wells]
    final_destinations = elution_plate.rows()[0][:Elution.num_wells]

    #Lysis.reagent_reservoir = reagent_res.rows()[0][:Lysis.num_wells] # 1 row, 4 columns (first ones)
    #Beads_PK.reagent_reservoir = reagent_res.rows()[0][Lysis.num_wells:(Lysis.num_wells+Beads_PK.num_wells)]
    #VHB.reagent_reservoir = reagent_res.rows()[0][(Lysis.num_wells+Beads_PK.num_wells):(Lysis.num_wells + Beads_PK.num_wells + VHB.num_wells)]
    #SPR.reagent_reservoir = reagent_res.rows()[0][VHB.num_wells:(Lysis.num_wells + Beads_PK.num_wells + VHB.num_wells + SPR.num_wells)]
    #SPR.reagent_reservoir = reagent_res_2.rows()[0][:SPR.num_wells]
    #Water.reagent_reservoir = reagent_res_2.rows()[0][-1]
    #work_destinations = deepwell_plate.rows()[0][:Elution.num_wells]
    #final_destinations = elution_plate.rows()[0][:Elution.num_wells]

    # pipettes. P1000 currently deactivated
    m300 = ctx.load_instrument('p300_multi_gen2', 'left', tip_racks=tips300) # Load multi pipette
    #p1000 = ctx.load_instrument('p1000_single_gen2', 'left', tip_racks=tips1000) # load P1000 pipette

    #### used tip counter and set maximum tips available
    tip_track = {
        'counts': {m300: 0},
        'maxes': {m300: 96 * len(m300.tip_racks)} #96 tips per tiprack * number or tipracks in the layout
        }
        #, p1000: len(tips1000)*96}

###############################################################################

    ###############################################################################
    # STEP 1 MIX BEADS
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
    ### PREMIX BEADS
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        if not m300.hw_pipette['has_tip']:
            pick_up(m300) #These tips are reused in the first transfer of beads
            ctx.comment(' ')
            ctx.comment('Tip picked up')
        ctx.comment(' ')
        ctx.comment('Mixing '+ Beads.name)
        ctx.comment(' ')
        #Mixing
        custom_mix(m300, Beads, Beads.reagent_reservoir[Beads.col], vol = 180,
        rounds = 20, blow_out = False, mix_height = 0, offset = 0)
        ctx.comment('Finished premixing!')
        ctx.comment('Now, reagents will be transferred to deepwell plate.')
        ctx.comment(' ')
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))

    ###############################################################################
    # STEP 2 TRANSFER LYSIS
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
    #Transfer lysis
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        # max_volume_allowed = 160 # Tips allow up to 200uL, but we only allow max_volume_allowed
        lysis_trips = math.ceil(Lysis.reagent_volume / Lysis.max_volume_allowed)
        lysis_volume = Lysis.reagent_volume / lysis_trips #136.66
        lysis_transfer_vol = []
        for i in range(lysis_trips):
            lysis_transfer_vol.append(lysis_volume + Lysis.disposal_volume)
        #ctx.comment(print(lysis_transfer_vol))
        #lysis_transfer_vol = [lysis_volume + Lysis.disposal_volume, lysis_volume + Lysis.disposal_volume, lysis_volume + Lysis.disposal_volume] # Three rounds of 140 + disposal volume
        x_offset_source = 0
        x_offset_dest   = 0
        rinse = True
        for i in range(num_cols):
            ctx.comment("Column: " + str(i))
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for j,transfer_vol in enumerate(lysis_transfer_vol):
                #Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(Lysis, multi_well_rack_area, transfer_vol * 8)

                if change_col == True: #If we switch column because there is not enough volume left in current reservoir column we mix new column
                    ctx.comment('Mixing new reservoir column: ' + str(Lysis.col))
                    custom_mix(m300, Lysis, Lysis.reagent_reservoir[Lysis.col],
                    vol = 180, rounds = 10, blow_out = False, mix_height = 0, offset = 0)
                ctx.comment('Aspirate from reservoir column: ' + str(Lysis.col))
                ctx.comment('Pickup height is ' + str(pickup_height))
                #if j!=0:
                #    rinse = False
                move_vol_multi(m300, reagent = Lysis, source = Lysis.reagent_reservoir[Lysis.col],
                dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 2, blow_out = False)
                #m300.move_to(work_destinations[i].top(0))
                #m300.air_gap(Lysis.air_gap_vol_bottom) #air gap
            #m300.aspirate(air_gap_vol_bottom, work_destinations[i].top(10), rate = Lysis.flow_rate_aspirate) #air gap
            #m300.dispense(disposal_volume + air_gap_vol_bottom, location = Lysis.reagent_reservoir[Lysis.col].top(0), rate = Lysis.flow_rate_dispense)
            ctx.comment(' ')
            ctx.comment('Mixing sample ')
            custom_mix(m300, Lysis, location = work_destinations[i], vol = 180,
            rounds = 20, blow_out = False, mix_height = 0, offset = 0)
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(Lysis.air_gap_vol_bottom) #air gap
            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))

    ###############################################################################
    # STEP 3 INCUBATING WITHOUT MAGNET
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
    #Transfer magnetic beads
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating for ' + format(STEPS[STEP]['wait_time']) + ' seconds.') # minutes=2
        ctx.comment(' ')
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 3 TRANSFER MAGNET BEADS
        ########

    ###############################################################################
    # STEP 4 INCUBATE WAIT WITH MAGNET ON
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        magdeck.engage(height=mag_height)
        ctx.comment(' ')
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating ON magnet for ' + format(STEPS[STEP]['wait_time']) + ' seconds.') # minutes=2
        ctx.comment(' ')
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 4 INCUBATE WAIT WITH MAGNET ON
        ########

    ###############################################################################
    # STEP 5 REMOVE SUPERNATANT
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        # remove supernatant -> height calculation can be omitted and referred to bottom!
        #max_volume_allowed = 190 # Tips allow up to 200uL, but we only allow max_volume_allowed
        supernatant_trips = math.ceil((Lysis.reagent_volume + sample_volume) / Lysis.max_volume_allowed)
        supernatant_volume = Lysis.max_volume_allowed # We try to remove an exceeding amount of supernatant to make sure it is empty
        #supernatant_volume = (Lysis.reagent_volume + sample_volume) / supernatant_trips #136.66
        supernatant_transfer_vol = []
        for i in range(supernatant_trips):
            supernatant_transfer_vol.append(supernatant_volume + Elution.disposal_volume)
        #ctx.comment(print(supernatant_transfer_vol))
        #supernatant_volume = 150
        #supernatant_vol = [supernatant_volume + Elution.disposal_volume, supernatant_volume + Elution.disposal_volume, supernatant_volume + Elution.disposal_volume, supernatant_volume + Elution.disposal_volume]
        x_offset_rs = 2

        for i in range(num_cols):
            x_offset_source = find_side(i) * x_offset_rs
            x_offset_dest   = 0
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in supernatant_transfer_vol:
                #Pickup_height is fixed here
                pickup_height = 1 # Original 0.5
                ctx.comment('Aspirate from deep well column: ' + str(i+1))
                ctx.comment('Pickup height is ' + str(pickup_height) +' (fixed)')
                move_vol_multi(m300, reagent = Elution, source = work_destinations[i],
                dest = waste, vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 2, blow_out = False)
                #m300.air_gap(Elution.air_gap_vol_bottom) #air gap
            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 5 REMOVE SUPERNATANT
        ########

    ###############################################################################
    # STEP 6 MAGNET OFF
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch off magnet
        magdeck.disengage()

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 6 MAGNET OFF
        ########

    ###############################################################################
    # STEP 7 VHB/WB1
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        #max_volume_allowed = 190 # Tips allow up to 200uL, but we only allow max_volume_allowed
        vhb_trips = math.ceil(VHB.reagent_volume / VHB.max_volume_allowed)
        vhb_volume = VHB.reagent_volume / vhb_trips #136.66
        vhb_transfer_vol = []
        for i in range(vhb_trips):
            vhb_transfer_vol.append(vhb_volume + VHB.disposal_volume)
        #vhb_volume = 166.66
        #vhb_wash_vol = [vhb_volume + VHB.disposal_volume, vhb_volume + VHB.disposal_volume, vhb_volume + VHB.disposal_volume]
        x_offset_rs = 2.5
        rinse = False # Not needed

        ########
        # whb washes
        for i in range(num_cols):
            x_offset_source = 0
            x_offset_dest   = -1 * find_side(i) * x_offset_rs
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in vhb_transfer_vol:
                #Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(VHB, multi_well_rack_area, transfer_vol*8)
                ctx.comment('Aspirate from Reservoir column: ' + str(VHB.col))
                ctx.comment('Pickup height is ' + str(pickup_height))
                #if i!=0:
                #    rinse = False
                move_vol_multi(m300, reagent = VHB, source = VHB.reagent_reservoir[VHB.col],
                dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 2, blow_out = False)
                #m300.move_to(work_destinations[i].top(0))
                #m300.air_gap(VHB.air_gap_vol_bottom) #air gap
                #m300.aspirate(air_gap_vol_bottom, work_destinations[i].top(10), rate = VHB.flow_rate_aspirate) #air gap
            custom_mix(m300, VHB, location = work_destinations[i], vol = 180,
                rounds = 20, blow_out = False, mix_height = 0, offset = x_offset_dest)
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(VHB.air_gap_vol_bottom) #air gap
            #m300.aspirate(air_gap_vol_bottom, work_destinations[i].top(10), rate = VHB.flow_rate_aspirate) #air gap
            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 7 ADD VHB
        ########

    ###############################################################################
    # STEP 8 INCUBATE WAIT WITH MAGNET ON
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        # switch on magnet
        magdeck.engage(mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Wait for 5 minutes.')
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ####################################################################
        # STEP 8* WAIT FOR 5'
        ########

    ###############################################################################
    # STEP 9 REMOVE SUPERNATANT
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        # remove supernatant -> height calculation can be omitted and referred to bottom!
        #max_volume_allowed = 190 # Tips allow up to 200uL, but we only allow max_volume_allowed
        supernatant_trips = math.ceil(VHB.reagent_volume / VHB.max_volume_allowed)
        supernatant_volume = VHB.max_volume_allowed # We try to remove an exceeding amount of supernatant to make sure it is empty
        #supernatant_volume = VHB.reagent_volume / supernatant_trips #136.66
        supernatant_transfer_vol = []
        for i in range(supernatant_trips):
            supernatant_transfer_vol.append(supernatant_volume + Elution.disposal_volume)
        #supernatant_vol = [175, 175, 175, 175, 175, 175, 175]
        x_offset_rs = 2

        for i in range(num_cols):
            x_offset_source = find_side(i) * x_offset_rs
            x_offset_dest   = 0
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in supernatant_transfer_vol:
                #Pickup_height is fixed here
                pickup_height = 1 # Original 0.5
                ctx.comment('Aspirate from deep well column: ' + str(i+1))
                ctx.comment('Pickup height is ' + str(pickup_height) +' (fixed)')
                move_vol_multi(m300, reagent = Elution, source = work_destinations[i],
                dest = waste, vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 2, blow_out = False)
                #m300.air_gap(Elution.air_gap_vol_bottom) #air gap
            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 9 REMOVE SUPERNATANT
        ########

    ###############################################################################
    # STEP 10 MAGNET OFF
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch off magnet
        magdeck.disengage()

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+ str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 10 MAGNET OFF
        ########

    ###############################################################################
    # STEP 11 SPR
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        #max_volume_allowed = 190 # Tips allow up to 200uL, but we only allow max_volume_allowed
        spr_trips = math.ceil(SPR.reagent_volume / SPR.max_volume_allowed)
        spr_volume = SPR.reagent_volume / spr_trips #136.66
        spr_transfer_vol = []
        for i in range(spr_trips):
            spr_transfer_vol.append(spr_volume + SPR.disposal_volume)
        #spr_volume = 166.66
        #spr_wash_vol = [spr_volume + SPR.disposal_volume, spr_volume + SPR.disposal_volume, spr_volume + SPR.disposal_volume]
        x_offset_rs = 2.5
        rinse = False # No rinse needed

        ########
        # spr washes
        for i in range(num_cols):
            x_offset_source = 0
            x_offset_dest   = -1 * find_side(i) * x_offset_rs
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in spr_transfer_vol:
                #Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(SPR, multi_well_rack_area, transfer_vol*8)
                ctx.comment('Aspirate from Reservoir column: ' + str(SPR.col))
                ctx.comment('Pickup height is ' + str(pickup_height))
                #if i!=0:
                #    rinse = False
                move_vol_multi(m300, reagent = SPR, source = SPR.reagent_reservoir[SPR.col],
                dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 2, blow_out = False)
                #m300.move_to(work_destinations[i].top(0))
                #m300.air_gap(SPR.air_gap_vol_bottom) #air gap
                #m300.aspirate(air_gap_vol_bottom, work_destinations[i].top(10), rate = SPR.flow_rate_aspirate) #air gap
            custom_mix(m300, VHB, location = work_destinations[i], vol = 180,
                rounds = 20, blow_out = False, mix_height = 0, offset = x_offset_dest)
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(SPR.air_gap_vol_bottom) #air gap
            #m300.aspirate(air_gap_vol_bottom, work_destinations[i].top(10), rate = SPR.flow_rate_aspirate) #air gap
            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 11 ADD SPR
        ########

    ###############################################################################
    # STEP 12 INCUBATE WAIT WITH MAGNET ON
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        # switch on magnet
        magdeck.engage(mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Wait for 5 minutes.')
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ####################################################################
        # STEP 12* WAIT FOR 5'
        ########

    ###############################################################################
    # STEP 13 REMOVE SUPERNATANT
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        # remove supernatant -> height calculation can be omitted and referred to bottom!
        #max_volume_allowed = 190 # Tips allow up to 200uL, but we only allow max_volume_allowed
        supernatant_trips = math.ceil(SPR.reagent_volume / SPR.max_volume_allowed)
        supernatant_volume = SPR.max_volume_allowed # We try to remove an exceeding amount of supernatant to make sure it is empty
        #supernatant_volume = SPR.reagent_volume / supernatant_trips #136.66
        supernatant_transfer_vol = []
        for i in range(supernatant_trips):
            supernatant_transfer_vol.append(supernatant_volume + Elution.disposal_volume)
        #supernatant_vol = [175, 175, 175, 175, 175]
        x_offset_rs = 2

        for i in range(num_cols):
            x_offset_source = find_side(i) * x_offset_rs
            x_offset_dest   = 0
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in supernatant_transfer_vol:
                #Pickup_height is fixed here
                pickup_height = 1 # Original 0.5
                ctx.comment('Aspirate from deep well column: ' + str(i+1))
                ctx.comment('Pickup height is ' + str(pickup_height) +' (fixed)')
                move_vol_multi(m300, reagent = Elution, source = work_destinations[i],
                dest = waste, vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 2, blow_out = False)
                #m300.air_gap(Elution.air_gap_vol_bottom) #air gap
            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 13 REMOVE SUPERNATANT
        ########

    ###############################################################################
    # STEP 14 MAGNET OFF
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch off magnet
        magdeck.disengage()

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 14 MAGNET OFF
        ########

    ###############################################################################
    # STEP 15 SPR
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        #max_volume_allowed = 190
        spr_trips = math.ceil(SPR.reagent_volume / SPR.max_volume_allowed)
        spr_volume = SPR.reagent_volume / spr_trips #136.66
        spr_transfer_vol = []
        for i in range(spr_trips):
            spr_transfer_vol.append(spr_volume + SPR.disposal_volume)
        #spr_volume = 166.66
        #spr_wash_vol = [spr_volume + SPR.disposal_volume, spr_volume + SPR.disposal_volume, spr_volume + SPR.disposal_volume]
        x_offset_rs = 2.5
        rinse = False #No rinse needed
        ########
        # spr washes
        for i in range(num_cols):
            x_offset_source = 0
            x_offset_dest   = -1 * find_side(i) * x_offset_rs
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in spr_transfer_vol:
                #Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(SPR, multi_well_rack_area, transfer_vol*8)
                ctx.comment('Aspirate from Reservoir column: ' + str(SPR.col))
                ctx.comment('Pickup height is ' + str(pickup_height))
                #if i!=0:
                #    rinse = False
                move_vol_multi(m300, reagent = SPR, source = SPR.reagent_reservoir[SPR.col],
                dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                pickup_height = pickup_height, rinse = rinse, avoid_droplet = False, wait_time = 2, blow_out = False)
                #m300.aspirate(SPR.air_gap_vol_bottom, work_destinations[i].top(10), rate = SPR.flow_rate_aspirate) #air gap
            custom_mix(m300, VHB, location = work_destinations[i], vol = 180,
                rounds = 20, blow_out = False, mix_height = 0, offset = x_offset_dest)
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(SPR.air_gap_vol_bottom) #air gap
            #m300.aspirate(air_gap_vol_bottom, work_destinations[i].top(10), rate = SPR.flow_rate_aspirate) #air gap
            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 15 ADD SPR
        ########

    ###############################################################################
    # STEP 16 INCUBATE WAIT WITH MAGNET ON
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        # switch on magnet
        magdeck.engage(mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Wait for 30 seconds.')
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ####################################################################
        # STEP 16* WAIT FOR 5'
        ########

    ###############################################################################
    # STEP 17 REMOVE SUPERNATANT
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        # remove supernatant -> height calculation can be omitted and referred to bottom!
        #max_volume_allowed = 190 # Tips allow up to 200uL, but we only allow max_volume_allowed
        supernatant_trips = math.ceil(SPR.reagent_volume / SPR.max_volume_allowed)
        supernatant_volume = SPR.max_volume_allowed # We try to remove an exceeding amount of supernatant to make sure it is empty
        #supernatant_volume = SPR.reagent_volume / supernatant_trips #136.66
        supernatant_transfer_vol = []
        for i in range(supernatant_trips):
            supernatant_transfer_vol.append(supernatant_volume + Elution.disposal_volume)
        #supernatant_vol = [175, 175, 175]
        x_offset_rs = 2

        for i in range(num_cols):
            x_offset_source = find_side(i) * x_offset_rs
            x_offset_dest   = 0
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in supernatant_transfer_vol:
                #Pickup_height is fixed here
                pickup_height = 1 # Original 0.5
                ctx.comment('Aspirate from deep well column: ' + str(i+1))
                ctx.comment('Pickup height is ' + str(pickup_height) +' (fixed)')
                move_vol_multi(m300, reagent = Elution, source = work_destinations[i],
                dest = waste, vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 2, blow_out = False)
                #m300.move_to(waste.top(0))
                #m300.air_gap(Elution.air_gap_vol_bottom) #air gap
            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)
        ctx.comment('Used tips in total: ' + str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 17 REMOVE SUPERNATANT
        ########

    ###############################################################################
    # STEP 18 ALLOW DRY
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Incubating OFF magnet for ' + format(STEPS[STEP]['wait_time']) + ' seconds.') # minutes=2
        ctx.comment(' ')
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)
        ctx.comment('Used tips in total: ' + str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 18 ALLOW DRY
        ########

    ###############################################################################
    # STEP 19 MAGNET OFF
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')

        # switch off magnet
        magdeck.disengage()

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)
        ctx.comment('Used tips in total: ' + str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 19 MAGNET OFF
        ########

    ###############################################################################
    # STEP 20 Transfer water
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        #Water elution
        #max_volume_allowed = 190 # Tips allow up to 200uL, but we only allow max_volume_allowed
        water_trips = math.ceil(Water.reagent_volume / Water.max_volume_allowed)
        water_volume = Water.reagent_volume / water_trips
        water_wash_vol = []
        for i in range(water_trips):
            water_wash_vol.append(water_volume + Elution.disposal_volume)
        #water_wash_vol = [50 + Water.disposal_volume]
        x_offset_rs = 2.5

        ########
        # Water or elution buffer
        for i in range(num_cols):
            x_offset_source = 0
            x_offset_dest   = -1 * find_side(i) * x_offset_rs # Original 0
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in water_wash_vol:
                #Calculate pickup_height based on remaining volume and shape of container
                [pickup_height, change_col] = calc_height(Water, multi_well_rack_area, transfer_vol*8)
                ctx.comment('Aspirate from Reservoir column: ' + str(Water.col))
                ctx.comment('Pickup height is ' + str(pickup_height))

                move_vol_multi(m300, reagent = Water, source = Water.reagent_reservoir,
                dest = work_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 0, blow_out = False)

                #m300.move_to(work_destinations[i].top(0))
                #m300.air_gap(Water.air_gap_vol_bottom) #air gap
                #m300.aspirate(air_gap_vol_bottom, work_destinations[i].top(10), rate = Water.flow_rate_aspirate) #air gap
            ctx.comment(' ')
            ctx.comment('Mixing sample with Water')
            #Mixing
            custom_mix(m300, Elution, work_destinations[i], vol = 40, rounds = 20,
            blow_out = False, mix_height = 0, offset = x_offset_dest)
            m300.move_to(work_destinations[i].top(0))
            m300.air_gap(Water.air_gap_vol_bottom) #air gap
            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
            tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 20 Transfer water
        ########

    ###############################################################################
    # STEP 21 WAIT FOR 10'
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Wait for ' + format(STEPS[STEP]['wait_time']) + ' seconds.')
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ####################################################################
        # STEP 21* WAIT FOR 10'
        ########

    ###############################################################################
    # STEP 22 INCUBATE WAIT WITH MAGNET ON
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        # switch on magnet
        magdeck.engage(mag_height)
        ctx.delay(seconds=STEPS[STEP]['wait_time'], msg='Wait for 5 minutes.')
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ####################################################################
        # STEP 22* WAIT FOR 5'
        ########

    ###############################################################################
    # STEP 23 TRANSFER TO ELUTION PLATE
    ########
    STEP += 1
    if STEPS[STEP]['Execute']==True:
        start = datetime.now()
        ctx.comment(' ')
        ctx.comment('###############################################')
        ctx.comment('Step '+str(STEP)+': '+STEPS[STEP]['description'])
        ctx.comment('###############################################')
        ctx.comment(' ')
        #max_volume_allowed = 150 # Tips allow up to 200uL, but we only allow max_volume_allowed
        elution_trips = math.ceil(Elution.reagent_volume / Elution.max_volume_allowed)
        elution_volume = Elution.reagent_volume / elution_trips
        elution_vol = []
        for i in range(elution_trips):
            elution_vol.append(elution_volume + Elution.disposal_volume)
        #elution_vol =[50]
        x_offset_rs = 2
        for i in range(num_cols):
            x_offset_source = find_side(i) * x_offset_rs
            x_offset_dest   = (find_side(i) * x_offset_rs)
            if not m300.hw_pipette['has_tip']:
                pick_up(m300)
            for transfer_vol in elution_vol:
                #Pickup_height is fixed here
                pickup_height = 1
                ctx.comment('Aspirate from deep well column: ' + str(i+1))
                ctx.comment('Pickup height is ' + str(pickup_height) +' (fixed)')

                move_vol_multi(m300, reagent = Elution, source = work_destinations[i],
                dest = final_destinations[i], vol = transfer_vol, x_offset_source = x_offset_source, x_offset_dest = x_offset_dest,
                pickup_height = pickup_height, rinse = False, avoid_droplet = False, wait_time = 2, blow_out = False)

                #m300.move_to(final_destinations[i].top(0))
                #m300.air_gap(Elution.air_gap_vol_bottom) #air gap
            if recycle_tip == True:
                m300.return_tip()
            else:
                m300.drop_tip(home_after = False)
                tip_track['counts'][m300] += 8
        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('Step ' + str(STEP) + ': ' + STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:']=str(time_taken)
        ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
        ###############################################################################
        # STEP 23 TRANSFER TO ELUTION PLATE
        ########

    '''if not ctx.is_simulating():
        with open(file_path,'w') as outfile:
            json.dump(STEPS, outfile)'''

    ctx.comment(' ')
    ctx.comment('###############################################')
    ctx.comment('Homing robot')
    ctx.comment('###############################################')
    ctx.comment(' ')
    ctx.home()

    # Disengage magnet
    magdeck.disengage()
###############################################################################
    # Light flash end of program
    import os
    #os.system('mpg123 /etc/audio/speaker-test.mp3')
    for i in range(3):
        gpio.set_rail_lights(False)
        gpio.set_button_light(1,0,0)
        time.sleep(0.3)
        gpio.set_rail_lights(True)
        gpio.set_button_light(0,0,1)
        time.sleep(0.3)
    gpio.set_button_light(0,1,0)
    ctx.comment('Finished! \nMove deepwell plate (slot 5) to Station C for MMIX addition and PCR preparation.')
    ctx.comment('Used tips in total: '+str(tip_track['counts'][m300]))
    ctx.comment('Used racks in total: '+str(tip_track['counts'][m300]/96))
    ctx.comment('Available tips: '+str(tip_track['maxes'][m300]))
