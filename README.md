# Código para el simulador 2D de elementos finitos en Julia

En este repositorio se muestra el código creado para construir un simulador 2D de elementos finitos en Julia. A continuación se documentará la información más importante a tener en cuenta sobre el código.

# Tabla de Contenidos

- [Instalación](#instalación)
- [Uso](#uso)
- [Licencia](#licencia)

# Instalación
Para poder utilizar el código debes instalar varios paquetes de Julia e importar la  API de Gmsh. Los paquetes o librerias a instalar son: 
- LinearAlgebra: Para obtener los autovalores y autovectores de la matriz de masa y rigidez
- SimplexQuad (https://github.com/eschnett/SimplexQuad.jl)
- Plots: para representar las mallas, asi como los resultados finales
- LaTeXStrings: para personalizar de manera matemática nuestros gráficos

## API de Gmsh
Para utilizar la API de Gmsh necesitamos la biblioteca dinámica
Gmsh y el módulo Julia (‘gmsh.jl’). Esto se obtiene del kit binario de desarrollo de software (SDK) disponible en el sitio web de Gmsh (https://gmsh.info/bin/), para Windows, Linux y macOS, descargando y descomprimiendo el archivo gmsh*-sdk.* correspondiente a su sistema operativo y agregando el directorio
"lib"del SDK a JULIA_LOAD_PATH.  Una vez tengamos descargado y descomprimido el SDK, buscamos la carpeta /.julia/config,
en caso de que no exista, tendremos que crearla. Una vez esté la carpeta creada, hay
que crear o editar el archivo startup.jl de la siguiente forma:

```Julia
 # Guardamos la ruta al archivo que queremos cargar en
JULIA_LOAD_PAD push !( LOAD_PATH , "/ ruta/al/archivo ")
```

Este archivo startup.jl se utiliza para realizar configuraciones personalizadas y cargar paquetes o módulos específicos que deseas que estén disponibles cada vez que
inicies Julia. Como último paso habría que reiniciar Julia.

# Uso
El código consta de una serie de modulos, donde el principal es el denominado "Numerical_integration_v1.jl". Este se encarga de llamar a "Mesh_object.jl" que a su vez llama a todos los demás. Las dependencias de los modulos se muestran en la siguiente imagen:

![dependencias_código](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/2b5e0d2e-07d6-4e69-8abd-a097c1722f26)

Donde las flechas muestran de quién obtine la información o de quién depende.

## Rect_waveguide_triangles/quads.jl 
Estas dos fracciones de código utilizan la API de Gmsh en Julia para la creación del objeto a utilizar. De ellas obtenos información importante como los nodos, las coordenadas de los mismos, los grupos físicos (condiciones de contorno), etc.

## Elements.jl
Este módulo se utiliza para guardar los diferentes tipos de elementos. Se han implementado ambos elementos, triangulos y quads. Aun asi no se ha realizado l orientación de quads, que queda pendiente para futuros trabajos. En este modulo se llama a "Elements_order.jl" donde se guarda el orden de los elementos transformandolo a un formato legible a partir de los códigos internos de Gmsh mostrados en esta página web: https://docs.juliahub.com/GmshTools/9rYp5/0.4.2/element_types/ .

## Material.jl
Gracias a este módulo se pueden guardar los diversos materiales de los que se compone la malla. Será util para futuros trabajos donde se analicen mallas genéricas. De momento no tiene influencia en el código.

## Boundary_cond.jl
En este módulo guardamos las condiciones de contorno presentes en la malla asi como todos los nodos que la sufren. Este modulo se utiliza para poder comprobar que nodos del elemento tienen determinada condición de contorno.

## Mesh_object.jl
