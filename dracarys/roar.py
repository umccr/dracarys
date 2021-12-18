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
      \ `\#    #/  -.'`_///`_///
       ;  \#  #/ '.__.'=\_.'
       |   '-#;    _|====\_
       ;      '  /`  `\==| |
        \     .        \=| /
         '-.._         // /__
              `)-.    `----._|
              <_________\_\_|
"""
