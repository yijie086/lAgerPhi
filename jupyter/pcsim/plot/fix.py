## various fixes for ROOT figures
import ROOT

def fix_palette(histo):
    ## dummy draw to initiate the palette
    ROOT.gPad.SetRightMargin(.15)
    histo.draw('colz')
    palette = histo.GetListOfFunctions().FindObject('palette')
    palette.SetX1NDC(.86)
    palette.SetX2NDC(.9)
    ROOT.gPad.Draw()
    ROOT.gPad.Clear()
