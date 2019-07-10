# -*- coding: utf-8 -*-
"""
Created on Mon May 28 22:42:59 2018

@author: shamm
"""

def DSSStartup(): 
    import win32com.client
    try:
        DSSObj = win32com.client.Dispatch("OpenDSSEngine.DSS")
        DSSText = DSSObj.Text
        DSSCircuit = DSSObj.ActiveCircuit
        DSSSolution = DSSCircuit.Solution
        DSSElem = DSSCircuit.ActiveCktElement
        DSSBus = DSSCircuit.ActiveBus
        DSSText.command = 'clear'
    except Exception as inst:
        print(type(inst))
        return {}
    return {'dssobj': DSSObj, 'dsstext': DSSText, 'dsscircuit': DSSCircuit, 'dsssolution' : DSSSolution, 'dsselem' : DSSElem, 'dssbus': DSSBus   }