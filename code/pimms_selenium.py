from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager

import numpy as np
from astropy.table import Table

'''
options = webdriver.ChromeOptions()
options.add_argument('headless')
options.add_argument('window-size=1920x1080')
options.add_argument("disable-gpu")

driver = webdriver.Chrome(executable_path='/Users/minghao/Research/Softwares/chromedriver')
'''

def clean_and_update(ele, val):
    ele.clear()
    ele.send_keys(val)

def get_convfactor(redshift, logNH, Elow, Ehigh, cycle=3, mwnh=1e20):
    service = Service(ChromeDriverManager().install())
    options = webdriver.ChromeOptions()
    options.add_argument('headless')
#    options.add_argument('window-size=10*10')
    options.add_argument("disable-gpu")
    options.headless=True

    driver = webdriver.Chrome(service=service, options=options)

    # Navigate to the webpage
    driver.get("https://cxc.harvard.edu/toolkit/pimms.jsp")

    # Select the 'Count Rate' radio button under 'Input'
    driver.find_element(By.XPATH, "//*[@id=\"jspForm\"]/center[1]/table[1]/tbody/tr/td[1]/table/tbody/tr[1]/td/table/tbody/tr/td[1]/input").click()
    driver.find_element(By.XPATH, "//*[@id=\"jspForm\"]/center[1]/table[1]/tbody/tr/td[2]/table/tbody/tr[1]/td/table/tbody/tr/td[2]/input").click()

    # Select a mission from the dropdown
    Select(driver.find_element(By.NAME, "inputMissionSelector")).select_by_visible_text("CHANDRA-Cycle %d"%cycle)

    # Fill in the 'Input Energy' fields
    clean_and_update(driver.find_element(By.NAME, "inputEnergyLow"),"%.1f"%Elow)
    clean_and_update(driver.find_element(By.NAME, "inputEnergyHigh"),"%.1f"%Ehigh)
    clean_and_update(driver.find_element(By.NAME, "outputFluxEnergyLow"),"%.1f"%Elow)
    clean_and_update(driver.find_element(By.NAME, "outputFluxEnergyHigh"),"%.1f"%Ehigh)

    # Select the 'Power Law' model from the dropdown
    Select(driver.find_element(By.NAME, "modelSelector")).select_by_visible_text("Power Law")

    # Fill in the 'Galactic NH' field
    driver.find_element(By.NAME, "NH").send_keys("%s"%mwnh)
    driver.find_element(By.NAME, "redshift").send_keys("5")
    driver.find_element(By.NAME, "redshiftedNH").send_keys("1e%d"%logNH)
    driver.find_element(By.NAME, "photonIndex").send_keys("1.8")
    driver.find_element(By.NAME, "countRate").send_keys("1")


    # Click on the 'Calculate' button
    driver.find_element(By.ID, "calcButton").click()
    driver.implicitly_wait(10)
    # Assuming there is a result field with ID 'result'
    # Retrieve the output value
    output = driver.find_element(By.XPATH, "//*[@id=\"jspForm\"]/center[2]/table/tbody/tr/td/input").get_attribute('value')

    print("PIMMS Prediction:", output)

    driver.quit()
    return float(output)

def main():
    tbl = Table.read('/Users/minghao/Research/Projects/JWST/LRDs/data/new/all_lrds_cstack_good_mwnh_ff.fits')
    pimms_s_21 = []
    pimms_h_21 = []
    pimms_a_21 = []

    pimms_s_22 = []
    pimms_h_22 = []
    pimms_a_22 = []

    pimms_s_23 = []
    pimms_h_23 = []
    pimms_a_23 = []

    for index in range(len(tbl)):
        redshift = tbl['redshift'][index]
        name = tbl['names'][index]
        ra = tbl['RA'][index]
        dec = tbl['Dec'][index]
        if ra>200:#AEGIS
            cycle = 9
        elif ra>150:#CDFN
            cycle = 3
        else:# CDFS
            cycle = 10

        mwnh = tbl['NH_MW'][index]
#        print(str(mwnh))
#        continue
#        print(redshift, name, ra, dec, cycle)
#        '''
        pimms_s_21.append(get_convfactor(redshift, 21, 0.5, 2, cycle=cycle, mwnh=mwnh))
        pimms_h_21.append(get_convfactor(redshift, 21, 2, 8, cycle=cycle, mwnh=mwnh))
        pimms_a_21.append(get_convfactor(redshift, 21, 0.5, 8, cycle=cycle, mwnh=mwnh))

        pimms_s_22.append(get_convfactor(redshift, 22, 0.5, 2, cycle=cycle, mwnh=mwnh))
        pimms_h_22.append(get_convfactor(redshift, 22, 2, 8, cycle=cycle, mwnh=mwnh))
        pimms_a_22.append(get_convfactor(redshift, 22, 0.5, 8, cycle=cycle, mwnh=mwnh))

        pimms_s_23.append(get_convfactor(redshift, 23, 0.5, 2, cycle=cycle, mwnh=mwnh))
        pimms_h_23.append(get_convfactor(redshift, 23, 2, 8, cycle=cycle, mwnh=mwnh))
        pimms_a_23.append(get_convfactor(redshift, 23, 0.5, 8, cycle=cycle, mwnh=mwnh))

        #pimms_s_25.append(get_convfactor(redshift, 25, 0.5, 2, cycle=cycle))
        #pimms_h_25.append(get_convfactor(redshift, 25, 2, 8, cycle=cycle))
    '''
        pimms_s_243.append(get_convfactor(redshift, 24.3, 0.5, 2, cycle=cycle))
        pimms_h_243.append(get_convfactor(redshift, 24.3, 2, 8, cycle=cycle))
    '''
    tbl['pimms_s_21'] = pimms_s_21
    tbl['pimms_h_21'] = pimms_h_21
    tbl['pimms_a_21'] = pimms_a_21

    tbl['pimms_s_22'] = pimms_s_22
    tbl['pimms_h_22'] = pimms_h_22
    tbl['pimms_a_22'] = pimms_a_22

    tbl['pimms_s_23'] = pimms_s_23
    tbl['pimms_h_23'] = pimms_h_23
    tbl['pimms_a_23'] = pimms_a_23

#    tbl['pimms_s_25'] = pimms_s_25
#    tbl['pimms_h_25'] = pimms_h_25
    '''
    tbl['pimms_s_243'] = pimms_s_243
    tbl['pimms_h_243'] = pimms_h_243
    '''
    tbl.write('/Users/minghao/Research/Projects/JWST/LRDs/data/new/all_lrds_cstack_good_nH_pimms.fits', overwrite=True)
#    get_convfactor(5, 22, 2, 10, cycle=3)
if __name__=='__main__':
    main()
