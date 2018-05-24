#!/usr/local/bin/python2 -O
import ipywidgets as widgets
from scipy import sqrt
from IPython.display import display
from pcsim.util import classdict 
import json
import os
import subprocess
from tempfile import NamedTemporaryFile

PCSIM_PATH='../BUILD/program/pcsim-lp_gamma'
#PCSIM_PATH='../DEBUG3/program/pcsim-lp_gamma'

hbox_layout = widgets.Layout(display='flex',
                    flex_flow='row',
                    align_items='stretch',
                    border='none')
vbox_layout = widgets.Layout(display='flex', padding='10px 15px 5px 15px', width='45%')

def width(w):
    return widgets.Layout(width='%s' % w + '%')

def make_row(*args):
    return widgets.HBox([a['widget'] for a in args])

def title_widget(title):
    w = widgets.VBox([widgets.HTML('<h3>%s</h3>' % title)],
                     layout=vbox_layout)
    return classdict(widget=w)

def text_widget(title, caption, **args):
    title = widgets.HTML('<h4>%s</h4>' % title)
    caption= widgets.HTMLMath(caption, layout=width(70))
    text = widgets.Text(**args)
    w = widgets.VBox([title, caption, text], layout=vbox_layout)
    return  classdict(title=title, caption=caption, text=text, widget=w)

def progress_widget(title, caption, **args):
    title = widgets.HTML('<h4>%s</h4>' % title)
    caption= widgets.HTMLMath(caption, layout=width(70))
    text = widgets.Text(**args)
    w = widgets.VBox([title, caption, text], layout=vbox_layout)
    return  classdict(title=title, caption=caption, text=text, widget=w)

def menu_widget(title, caption, **args):
    title = widgets.HTML('<h4>%s</h4>' % title)
    caption= widgets.HTMLMath(caption, layout=width(70))
    menu = widgets.Dropdown(**args)
    w = widgets.VBox([title, caption, menu], layout=vbox_layout)
    return  classdict(title=title, caption=caption, menu=menu, widget=w)

def button_widget(title, caption, **args):
    title = widgets.HTML('<h4>%s</h4>' % title)
    caption= widgets.HTMLMath(caption, layout=width(70))
    button = widgets.Button(**args)
    w = widgets.VBox([title, caption, button], layout=vbox_layout)
    return  classdict(title=title, caption=caption, button=button, widget=w)

def int_slider_widget(title, caption, **args):
    title = widgets.HTML('<h4>%s</h4>' % title)
    caption= widgets.HTMLMath(caption, layout=width(70))
    range1 = widgets.IntSlider(readout=False, **args)
    range2 = widgets.BoundedIntText(**args)
    l = widgets.jslink((range1, 'value'), (range2, 'value'))
    w = widgets.VBox([title, caption, range1, range2], layout=vbox_layout)
    return classdict(title=title, caption=caption, slider=range1, widget=w)

def float_slider_widget(title, caption, **args):
    title = widgets.HTML('<h4>%s</h4>' % title)
    caption= widgets.HTMLMath(caption, layout=width(70))
    range1 = widgets.FloatSlider(readout=False, **args)
    range2 = widgets.BoundedFloatText(**args)
    l = widgets.jslink((range1, 'value'), (range2, 'value'))
    w = widgets.VBox([title, caption, range1, range2], layout=vbox_layout)
    return classdict(title=title, caption=caption, slider=range1, widget=w)

def float_range_widget(title, caption, **args):
    title = widgets.HTML('<h4>%s</h4>' % title)
    caption= widgets.HTMLMath(caption, layout=width(70))
    rmin = widgets.BoundedFloatText(description='min', value=args['min'], **args)
    rmax = widgets.BoundedFloatText(description='max', value=args['max'], **args)
    w = widgets.VBox([title, caption, rmin, rmax], layout=vbox_layout)
    return classdict(title=title, caption=caption, min=rmin, max=rmax, widget=w)

def optional_float_range_widget(title, caption, **args):
    raw = float_range_widget(title, caption, **args)
    raw.min.disabled = True
    raw.max.disabled = True
    raw.enable = widgets.Checkbox(value=False, description='Enable')
    raw.widget.children = [raw.title, raw.caption, raw.enable, raw.min, raw.max]
    def handler(change):
        raw.min.disabled = change['old']
        raw.max.disabled = change['old']
    raw.enable.observe(handler, 'value')
    return raw

def param_widget(title, caption, **args):
    title = widgets.HTML('<h4>%s</h4>' % title)
    caption= widgets.HTMLMath(caption, layout=width(70))
    param=classdict()
    col = [title, caption]
    for name in args['param']:
        param[name] = widgets.FloatText(**args['param'][name])
        col.append(param[name])
    w = widgets.VBox(col, layout=vbox_layout)
    return classdict(title=title, caption=caption, param=param, widget=w)
    

def empty_widget():
    w = widgets.VBox([], layout=vbox_layout)
    return classdict(widget=w)

def update_display(vbox, title, text):
    title = widgets.HTML('<h4>%s</h4>' % title)
    text = widgets.HTMLMath(text)
    vbox.children = [title, text]

def status_widget(title, **args):
    title = widgets.HTML('<h4>%s</h4>' % title)
    current = widgets.HTMLMath('')
    finished = widgets.HTMLMath('')
    w = widgets.VBox([title, current, finished], layout=vbox_layout)
    sw = classdict(widget=w, title=title, current=current, finished=finished, 
        fname='', fname_old = [])
    def handler(change):
        fname = change['fname']
        stat = ''
        if 'running' in change and change['running']:
            stat ='Generating:'
        else:
            stat += 'Ready to generate:'
        current.value = ''.join(['<b>%s</b><ul>' % stat, 
                            '<li>%s.root</li>' % fname,
                            '<li>%s.json</li>' % fname,
                            '<li>%s.log</li>' % fname,
                            '</ul>'])
        if 'fname_done' in change and change['fname_done'] not in sw.fname_old:
            f = '<li>%s.root</li>' % change['fname_done']
            if (f not in sw.fname_old):
                sw.fname_old.append(f)
        if len(sw.fname_old) > 0:
            text = '<b>Previously generated:</b><ul>'
            text += ''.join(sw.fname_old)
            text += '</ul>'
            finished.value = text
    sw.handler = handler
    return sw

##########################################################################################
### MC SETTINGS
##########################################################################################
mc = classdict()

mc.name = text_widget(
    'Name', 
    'This will be used as a tag for all output files.', 
    value='demo')
mc.nevents = int_slider_widget(
    'Number of events',
    'Total number of generated events.',
    min=100, max=1e7, value=1e5)
mc.tag = text_widget(
    'Tag', 
    'Additional tag to describe this MC.',
    value='none')
mc.run = int_slider_widget(
    'Run Number',
    'The run number is also used as random seed.',
    min=1, max=99999, value=1)

#rows = [make_row(mc.name, mc.tag), 
#        make_row(mc.run, mc.nevents)]
#mc.box = widgets.VBox(rows)

##########################################################################################
### BEAM SETTINGS
##########################################################################################
## titles
beam = classdict()
beam.lepton = classdict()
beam.target = classdict()

beam.lepton.title = title_widget('Lepton Beam')
beam.target.title = title_widget('Target Beam')

beam.lepton.energy = float_slider_widget(
    'Energy', 
    'Select the beam energy in GeV',
    min=.5, max=1000, value=5)
beam.target.energy = float_slider_widget(
    'Energy', 
    'Select the beam energy in GeV',
    min=.938272, max=1000, value=100)

## CM energy display
beam.CM = empty_widget()
def update_CM(change):
    Eb = beam.lepton.energy.slider.value
    Et = beam.target.energy.slider.value
    s = (Eb + Et)**2 - (sqrt(Eb**2 - .000511**2) - sqrt(Et**2 - .938272**2))**2
    update_display(
        beam.CM.widget, 
        'Total CM energy', 
        '<i>&radic;s</i> = %f GeV' % sqrt(s))
update_CM(0)
beam.lepton.energy.slider.observe(update_CM)
beam.target.energy.slider.observe(update_CM)

## lepton particle type
beam.lepton.type = menu_widget(
    'Type',
    'Select the particle type.',
    value='11', 
    options={'e': '11', 'e+': '-11'})

## target particle type
beam.target.type = menu_widget(
    'Type',
    'Select the particle type.',
    value='2212', 
    options={'p': '2212'})

rows = [
    make_row(beam.lepton.title, beam.target.title),
    make_row(beam.lepton.type, beam.target.type),
    make_row(beam.lepton.energy, beam.target.energy),
    make_row(beam.CM)
]

beam.box = widgets.VBox(rows)

##########################################################################################
### PHOTON SETTINGS
##########################################################################################

photon = classdict()
photon.vphoton = classdict()
photon.bs = classdict()
photon.bs.rl = classdict()

## generator type
photon.type = menu_widget(
    'Photon type',
    'Select the photon generator type.',
    options={
        'Virtual photon beam': 'vphoton', 
        'Bremsstrahlung beam': 'bremsstrahlung'},
    value='vphoton')

## BS model
photon.bs.model = menu_widget(
    'Bremsstrahlung model',
    'Select the bremsstrahlung model.',
    options={
        'Approximate': 'approx', 
        'Parameterization': 'param', 
        'Flat': 'flat'},
   value='flat')

## RL pickers
photon.bs.rl.param = menu_widget(
    'Radiator thickness',
    'Select the radiator thickness model.',
    options={
        '1% RL': .01, 
        '5% RL': .05, 
        '10% RL':.1},
   value=.1)
photon.bs.rl.flat = empty_widget()
photon.bs.rl.approx = float_slider_widget(
    'Radiator thickness',
    'Select the radiator thickness model [%].',
    min=1, max=300, value=1)

## y range
photon.vphoton.yrange = float_range_widget(
    '<i>y</i> range',
    'Select the <i>y</i> range.',
    min=.01, max=1)

## Q2 range
photon.vphoton.Q2range = float_range_widget(
    '<i>Q</i><sup>2</sup> range', 
    'Select the <i>Q</i><sup>2</sup> range.',
    min=.01, max=1e10)

## E range
photon.bs.Erange = float_range_widget(
    '<i>E</i> range',
    'Select the photon <i>E</i> range in GeV.',
    min=.001, max=1e9)

## W range
photon.vphoton.Wrange = optional_float_range_widget(
   '<i>W</i> range (optional)',
    'Click the check box to enable this cut.',
    min=.938272, max=1e9)


photon.box = widgets.VBox([])

def update_photon(change):
    rows = []
    if change['new'] == 'vphoton':
        rows = [
            make_row(photon.type, photon.vphoton.yrange),
            make_row(photon.vphoton.Q2range, photon.vphoton.Wrange)
        ]
    elif change['new'] == 'bremsstrahlung' or change['new'] in photon.bs.rl:
        rows = [
            make_row(photon.type, photon.bs.Erange),
            make_row(photon.bs.model,photon.bs.rl[photon.bs.model.menu.value])
        ]
    photon.box.children = rows

update_photon({'new':'vphoton'})

photon.type.menu.observe(update_photon, 'value')
photon.bs.model.menu.observe(update_photon, 'value')

##########################################################################################
### GENERATOR SETTINGS
##########################################################################################
gen = classdict()
gen.brodsky = classdict()
gen.pc = classdict()

## generator type
gen.model = menu_widget(
    'Process type',
    'Select the process type.',
    options={
        'Brodsky t-channel VM production': 'brodsky_2vmX', 
        'Charmed Pentaquark Production': 'gaussian_qpq'},
   value='brodsky_2vmX')

## Brodsky vm type
gen.brodsky.vmtype = menu_widget(
    'Vector meson type',
    'Select the vector meson type.',
    value=443,
    options={
        'J/psi': 443, 'Upsilon': 553})

## Brodsky recoil type
gen.brodsky.recoiltype = menu_widget(
    'Recoil type',
    'Select the recoil baryon type.',
    value=2212, 
    options={ 'proton': 2212, 'Delta+': 2214})

## Brodsky photo-production parameters
gen.brodsky.param = param_widget(
    'Photo-production cross section',
    'Usefull presets:' 
        '<ul>'
        '<li>threshold: xx, xx, xx</li>'
        '<li>HERA J/psi: xx, xx, xx</li>'
        '<li>HERA Upsilon: xx, xx, xx</li>'
        '</ul>',
    param={
        'b': {'description': 'b', 'value': 1.13},
        'c2g': {'description': '2-gluon', 'value': 6.499e3},
        'c3g': {'description': '3-gluon', 'value': 0}
    })

# R parameterization
gen.brodsky.R = param_widget(
    '<i>R</i> parameterization settings',
    '<i>R</i> parameterization settings',
    param={
        'factor': {'description': 'factor', 'value': 2.164},
        'power': {'description': 'power', 'value': 2.131}
    })
gen.pc.R = gen.brodsky.R

## Dipole parameterization
gen.brodsky.dipole = param_widget(
    '\'Dipole\' FF settings',
    'Necessary to relate photo-production to lepto-production',
    param={
        'power': {'description': 'power', 'value': 2.575}
    })
gen.pc.dipole = gen.brodsky.dipole
gen.empty = empty_widget()

gen.box = widgets.VBox([])

def update_gen(change):
    rows = []
    if change['new'] == 'brodsky_2vmX':
        rows = [
            make_row(gen.model, gen.brodsky.R),
            make_row(gen.brodsky.vmtype, gen.brodsky.dipole),
            make_row(gen.brodsky.recoiltype, gen.brodsky.param)]
    elif change['new'] == 'gaussian_qpq':
        rows = [
            make_row(gen.model, gen.pc.R),
            make_row(gen.empty, gen.pc.dipole)]
    gen.box.children = rows

update_gen({'new':'brodsky_2vmX'})

gen.model.menu.observe(update_gen, 'value')

##########################################################################################
### DETECTOR SETTINGS
##########################################################################################
detector = classdict()

## generator type
detector.setup = menu_widget(
    'Detector',
    'Select the detector.',
    options={
        'No detector' : '4pi',
        'JLEIC FastMC' : 'jleic'},
    value='4pi')

rows = [make_row(detector.setup)]
detector.box = widgets.VBox(rows)

##########################################################################################
### MAIN DASHBOARD
##########################################################################################
run = classdict()
run.button = button_widget(
    'Run',
    'Click here to run the MC',
    disabled=False,
    description='Run MC',
    button_style='success',
    tooltip='Click me',
    icon='check')

run.odir = text_widget(
    'Output directory',
    'Leave empty to use the current working directory.',
    value='')

run.status = status_widget('Status')

rows = [make_row(mc.name, mc.tag),
        make_row(mc.run, mc.nevents),
        make_row(run.odir),
        make_row(run.button, run.status)]
run.box = widgets.VBox(rows)

def fname():
    s = ''
    if mc.tag.text.value and mc.tag.text.value != 'none':
        s = '.'.join([
            'pcsim-lp_gamma', 
            mc.name.text.value,
            detector.setup.menu.value,
            mc.tag.text.value,
            'run%05i-%i' % (mc.run.slider.value, mc.nevents.slider.value)])
    else:
        s = '.'.join([
            'pcsim-lp_gamma', 
            mc.name.text.value,
            detector.setup.menu.value,
            'run%05i-%i' % (mc.run.slider.value, mc.nevents.slider.value)])
    if len(run.odir.text.value) == 0:
        return s
    elif run.odir.text.value[-1] != '/':
        return run.odir.text.value + '/' + s
    return run.odir.text.value + s

run.status.handler({'fname': fname()})


##########################################################################################
### MAIN DASHBOARD
##########################################################################################
def dashboard():
    title = widgets.HTMLMath('<h1>PCSIM $(l/\\gamma,p)$ Dashboard</h1>')

    accordion = widgets.Accordion(children=
        [i.box for i in [beam, photon, gen, detector, run]])
    accordion.set_title(0, 'Beam Settings')
    accordion.set_title(1, 'Photon Settings')
    accordion.set_title(2, 'Generator Settings')
    accordion.set_title(3, 'Detector Settings')
    accordion.set_title(4, 'Run PCSIM')
    display(title, accordion)

##########################################################################################
### COMPOSE THE MAIN JSON
##########################################################################################
def run_mc(dummy):
    my_fname = fname()
    run.button.button.disabled = True
    run.button.button.button_style='danger'
    run.button.button.description = 'Processing...'
    run.button.button.icon = 'hourglass'

    run.status.handler({'fname': my_fname, 'running': True})
    config = {
    'mc': {
        'type': 'pcsim-lp_gamma',
        'run': mc.run.slider.value,
        'events': mc.nevents.slider.value,
        'tag': mc.tag.text.value,
        'generator': {
            'type': mc.name.text.value,
            'beam': {
                'type': 'primary',
                'particle_type': beam.lepton.type.menu.value,
                'dir': [0, 0, 1],
                'energy': beam.lepton.energy.slider.value
            },
            'target': {
                'type': 'primary',
                'particle_type': beam.target.type.menu.value,
                'dir': [0, 0, -1],
                'energy': beam.target.energy.slider.value
            },
            'photon': {
                'type': photon.type.menu.value,
                'E_range': [photon.bs.Erange.min.value,
                            photon.bs.Erange.max.value]
            },
            'process_0': {
                'type': gen.model.menu.value,
                'vm_type': gen.brodsky.vmtype.menu.value,
                'recoil_type': gen.brodsky.recoiltype.menu.value
            }
        },
        'detector': {
            'type': detector.setup.menu.value
        }
    }}
    jsph = config['mc']['generator']['photon']
    if photon.type.menu.value == 'bremsstrahlung':
        jsph['model'] = photon.bs.model.menu.value
        if photon.bs.model.menu.value == 'param':
            jsph['rl'] = photon.bs.rl.param.menu.value
        elif photon.bs.model.menu.value == 'approx':
            jsph['rl'] = photon.bs.rl.approx.slider.value
    elif photon.type.menu.value == 'vphoton':
        jsph['y_range'] = [photon.vphoton.yrange.min.value, photon.vphoton.yrange.max.value]
        #if photon.vphoton.Q2range.enable.value:
        jsph['Q2_range'] = [photon.vphoton.Q2range.min.value,
            photon.vphoton.Q2range.max.value]
        if photon.vphoton.Wrange.enable.value:
            jsph['W_range'] = [photon.vphoton.Wrange.min.value,
                                photon.vphoton.Wrange.max.value]
        del jsph['E_range']
    ## gen
    jsgen = config['mc']['generator']['process_0']
    if gen.model.menu.value == 'brodsky_2vmX':
        jsgen['photo_b'] = gen.brodsky.param.param.b.value
        jsgen['photo_c2g'] = gen.brodsky.param.param.c2g.value
        jsgen['photo_c3g'] = gen.brodsky.param.param.c3g.value
        jsgen['R_vm_c'] = gen.brodsky.R.param.factor.value
        jsgen['R_vm_n'] = gen.brodsky.R.param.power.value
        jsgen['dipole_n'] = gen.brodsky.dipole.param.power.value
    print(config)

    fconf = open('config.json', 'w')
    json.dump(config, fconf, indent=4, sort_keys=True)
    fconf.write('')
    fconf.close()

    odir = '.'
    if len(run.odir.text.value):
        odir = run.odir.text.value

    #print ' '.join([PCSIM_PATH, '-cconfig.json', '-o %s' % odir])
    #p = subprocess.call([PCSIM_PATH, '-c config.json', '-o %s' % odir])
    os.system(' '.join([PCSIM_PATH, '-cconfig.json', '-o %s' % odir]))

    os.unlink(fconf.name)
    run.button.button.disabled = False
    run.button.button.button_style='success'
    run.button.button.description = 'Run MC'
    run.button.button.icon = 'check'
    run.status.handler({'fname': my_fname, 'fname_done': my_fname})
    
run.button.button.on_click(run_mc)

