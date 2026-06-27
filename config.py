import sys
import os
try:
    import lifelib
    assert(os.path.realpath(os.path.dirname(os.path.dirname(lifelib.__file__))) == os.path.realpath(os.path.dirname(__file__)))
except (ImportError, AssertionError):
    sys.stderr.write('''Please install lifelib in the current working directory by running:
                     git clone https://gitlab.com/apgoucher/lifelib
                     ''')
    raise SystemExit(1)
rule = sys.argv[1] if len(sys.argv) > 1 else 'b3s23'
rule = lifelib.genera.sanirule(rule)
print('Configuring lifelib for rule '+rule+'...')
lifelib.autocompile.set_rules(rule)
print('Successfully configured lifelib.')
