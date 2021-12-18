import click

avail_colours = [
    'black',
    'red',
    'green',
    'yellow',
    'blue',
    'magenta',
    'cyan',
    'white',
    'reset',
    'bright_black',
    'bright_red',
    'bright_green',
    'bright_yellow',
    'bright_blue',
    'bright_magenta',
    'bright_cyan',
    'bright_white'
    ]

@click.command()
@click.option('-c', '--color', default='green', type=click.Choice(avail_colours), help='Dragon colour')
def roar(color):
    """G'day"""
    click.echo(click.style(dragon, fg=color))

dragon = """
                     ,     _,
                    #\`-"-'/
                   #/   o (o
                  #/ \__   '._ 
  ,_#_#          #/  /=/`-. _")
   '-.`\#       #/  /=(_.. `-`.
      \ `\#    #/  -.'`_\\\`_\\\
       ;  \#  #/ '.__.'=\_.'
       |   '-#;    _|====\_ 
       ;      '  /`  `\==| \
        \     .        \=| /
         '-.._         // /__
              `)-.    `----._\
              <_________\_\_\
"""
