# Código para el simulador 2D de elementos finitos en Julia

En este repositorio se muestra el código creado para construir un simulador 2D de elementos finitos en Julia. A continuación se documentará la información más importante a tener en cuenta sobre el código.

# Tabla de Contenidos

- [Instalación](#instalación)
- [Uso](#uso)
- [Futuras ampliaciones](#futurasAmpliaciones)

# Instalación
Para poder utilizar el código debes instalar varios paquetes de Julia e importar la  API de Gmsh. Los paquetes o librerías a instalar son: 
- LinearAlgebra: Para obtener los autovalores y autovectores de la matriz de masa y rigidez
- SimplexQuad (https://github.com/eschnett/SimplexQuad.jl): nos permite obtener los puntos de integración y los pessos según la cuadratura de Gauss-Radau.
- Plots: para representar las mallas, así como los resultados finales
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
El código consta de una serie de módulos, donde el principal es el denominado "Numerical_integration_v1.jl". Este se encarga de llamar a "Mesh_object.jl" que a su vez llama a todos los demás. Las dependencias de los módulos se muestran en la siguiente imagen:

![dependencias_código](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/2b5e0d2e-07d6-4e69-8abd-a097c1722f26)

Donde las flechas muestran de quién obtine la información o de quién depende.

## Rect_waveguide_triangles/quads.jl 
Estos dos fragmentos de código utilizan la API de Gmsh en Julia para la creación del objeto a utilizar. De ellas obtenos información importante como los nodos, las coordenadas de los mismos, los grupos físicos (condiciones de contorno), etc.

## Elements.jl
Este módulo se utiliza para guardar los diferentes tipos de elementos. Se han implementado ambos elementos, triangulos y quads. Aun así, no se ha realizado la orientación de quads,la cual que queda pendiente para futuros trabajos. En este módulo se llama a "Elements_order.jl" donde se guarda el orden de los elementos transformándolo a un formato legible a partir de los códigos internos de Gmsh mostrados en esta página web: https://docs.juliahub.com/GmshTools/9rYp5/0.4.2/element_types/ .

## Material.jl
Gracias a este módulo se pueden guardar los diversos materiales de los que se compone la malla. Será útil para futuros trabajos donde se analicen mallas genéricas. De momento no tiene influencia en el código.

## Boundary_cond.jl
En este módulo guardamos las condiciones de contorno presentes en la malla así como todos los nodos que la sufren. Este módulo se utiliza para poder comprobar que nodos del elemento tienen determinada condición de contorno.

## Mesh_object.jl
Aqui se obtiene toda la información importante para ensamblar la malla. Se crean los elementos donde se tienen en cuenta las condiciones de contorno, se les asigna los grados de libertad, ect. Se guardan también los materiales de los que esta compuesta la malla así como todas las coordenadas de los nodos que la componen.

Toda esta información será necesaría para realizar los cálculos en "Numerical_integration_v1.jl".

### Mesh_plot.jl
En este módulo se crea un gráfico de la malla ensamblada utilizando toda la información creada y obtenida en "Mesh_object.jl". A continuación se muestran ejemplos de algunas mallas, de orden 1 y orden 2:

![triangular_mesh_order_1_boundary_1](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/4d66a5eb-a339-4931-8c25-b615d9f5d46b)

![triangular_mesh_order_2_boundary_1](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/2c6c1603-00d9-4976-86e7-a3ff2871cfb7)

## Numerical_integration_v1.jl
Este es el módulo main, ya que es donde se realizan todos los cálculos con la información obtenida hasta ahora y donde se obtienen los resultados finales. En el se calculan las matrices de masa y rigidez utilizadas para obtener los autovalores. Se puede elegir mediante el atributo 'FEM_TE' si se quiere calcular los modos TE o modos TM. 

Con los autovalores obtenidos, calculamos el error relativo, el cual disminuye a medida que aumentamos el orden de los elementos y hacemos la malla más fina con el atributo 'tm' del módulo "Rect_waveguide_triangles.jl". Debemos tener en cuenta dos puntos a la hora de calcular el error relativo en función de si calculamos modos TE o TM:
- Modos TE: el primer valor obtenido en 'k_c_fem' se utiliza como referencia, el segundo es el correspondiente al modo TE_10
  ![resultado_modos_TE](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/b73a9988-cc92-4e41-baee-d568cc5c5d65)

- Modos TM: no existen modos TM_10 o TM_01 por lo que el primer valor obtenido en 'k_c_fem' es el correspondiente al modo TM_11
  ![resultado_modos_TM](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/eaf11dc1-6092-47f8-9fae-81f0b001d839)

### Plot_results.jl
Este módulo se utiliza para graficar los resultados y realizar comprobaciones sobre ellos. Aqui se ha comprobado que la tasa de convergencia del error es O(h^2p).

![Mode_TE](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/159bb6d6-a2c8-4664-ac10-9a783ff31b60)

![Mode_TM](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/b3a7919a-edba-4cfd-aba0-0469e59ad6ed)

# Futuras Ampliaciones
El código se puede ampliar implementando las partes del código dedicadas a los cuadriláteros, así como implementando el uso de los materiales para poder resolver guias genéricas.
