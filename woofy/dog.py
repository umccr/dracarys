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
@click.option('-c', '--color', default='green', type=click.Choice(avail_colours), help='Dog colour')
def woofy(color):
    """G'day"""
    click.echo(click.style(dog, fg=color))

dog = """

                ______    ____
               :  ;;;;\__/:  ;|
                \;__/.... :  _/;
               ___:__ ..__\_/__;
               |  ## `--'   ##|;
               |_____/~;\_____|;
                 /~~~_ _ ~~   /
                 // (_:_)     |
           _     // ,'~ `,_||~||_
          //     ~~`,---,'~~~  |||
 ___     //         ~~~~      ;; \_  __
/_\/____::_        ,(:;:  __    ;; \/;;\  __
\_/) _  :: (       ; ;;:    \    / ;:;;::-,-'
   |[-]_::-|       : :;;:   /  * | ;:;;:;'
   | :__:: |       :.`,:::  : +  : /___:'
   |[_ ] [\|       ;. ;--`:_:.  *| /   /
   |__| |_]|    ,-' . :uu-'/     \|UUU/
   |_______|   ;_|_|_/    :_;_;_;_:
    [=====]

"""

