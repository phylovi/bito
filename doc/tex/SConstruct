import glob
import os

VariantDir('_build', src_dir='.')

env = Environment(ENV=os.environ)
inkscape = Builder(action = 'inkscape --without-gui --export-pdf=$TARGET $SOURCE')
env['BUILDERS']['Inkscape'] = inkscape
env['BUILDERS']['Latexdiff'] = Builder(action = 'latexdiff $SOURCES > $TARGET')
env['BUILDERS']['Copier'] = Builder(action = Copy('$TARGET', '$SOURCE'))

figure_pdfs = [env.Inkscape(target="figures/" + os.path.basename(svg).replace('.svg','.pdf'), source=svg)
               for svg in glob.glob('prep-figures/*.svg')]

pdfs = [env.Copier(target = '_build/' + os.path.basename(pdf), source = pdf)
        for pdf in glob.glob('figures/*.pdf')]

Depends(Flatten([pdfs]), Flatten([figure_pdfs]))

main=env.PDF(target='_build/main.pdf',source='main.tex')

Depends(Flatten([main]),
        Flatten([pdfs, figure_pdfs, 'main.bib']))

cont_build = env.Command('.continuous', ['main.bib', 'main.tex'],
    'while :; do inotifywait -e modify $SOURCES; scons -Q; done')
Alias('continuous', cont_build)

Default(main)
