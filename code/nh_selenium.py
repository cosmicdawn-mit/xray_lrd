from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager

import numpy as np
from astropy.table import Table

def clean_and_update(ele, val):
    ele.clear()
    ele.send_keys(val)

def get_nh_milkyway(ra, dec):
    service = Service(ChromeDriverManager().install())
    options = webdriver.ChromeOptions()
    options.add_argument('headless')
    options.add_argument('window-size=1920x1080')
    options.add_argument("disable-gpu")

    driver = webdriver.Chrome(service=service)

    # Navigate to the webpage
    driver.get("https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl")

    # input the coordinate
    clean_and_update(driver.find_element(By.XPATH, "/html/body/div[2]/div/div[5]/div/div/div[2]/table/tbody/tr/td[2]/center/form/table/tbody/tr[1]/td[2]/input"),"%.5f, %.5f"%(ra, dec))

    driver.find_element(By.XPATH, "/html/body/div[2]/div/div[5]/div/div/div[2]/table/tbody/tr/td[2]/center/form/p/input[1]").click()

#    clean_and_update(driver.find_element(By.ID, "/html/body/div[2]/div/div[5]/div/div/div[2]/table/tbody/tr/td[2]/center/form/table/tbody/tr[6]/td[2]/input"), "0.01")
    clean_and_update(driver.find_element(By.ID, "Radius"), "0.01")

    # Click on the 'Calculate' button
    driver.find_element(By.XPATH, "/html/body/div[2]/div/div[5]/div/div/div[2]/table/tbody/tr/td[2]/center/form/p/input[1]").click()
    driver.implicitly_wait(10)
    # Assuming there is a result field with ID 'result'
    # Retrieve the output value
    output = driver.find_element(By.XPATH, "/html/body/div[2]/div/div[5]/div/div/div[2]/table/tbody/tr/td[2]/center/pre/b[2]")

    textlist = output.text.split(' ')
    print("NH value:", float(textlist[-1]))

    driver.quit()
    return float(textlist[-1])
#    return float(output)

def main():
    # test the code
    tbl_all_lrds = Table.read('/Users/minghao/Research/Projects/JWST/LRDs/data/all_lrds_cstack_good.fits')
    mwhn_list = []
    for index in range(len(tbl_all_lrds)):
        ra = tbl_all_lrds['RA'][index]
        dec = tbl_all_lrds['Dec'][index]

        print(ra, dec)

        mwhn_list.append(get_nh_milkyway(ra, dec))

    tbl_all_lrds['NH_MW'] = mwhn_list

    tbl_all_lrds.write('/Users/minghao/Research/Projects/JWST/LRDs/data/new/all_lrds_cstack_good_mwnh.fits')



if __name__=='__main__':
    main()
