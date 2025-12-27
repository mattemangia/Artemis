# Artemis Physics Engine

Un motore fisico 2D leggero scritto in C# con demo interattiva.

## Struttura del Progetto

### ArtemisEngine
Motore fisico 2D che include:
- Sistema di corpi rigidi (RigidBody)
- Forme fisiche (cerchi e rettangoli)
- Rilevamento e risoluzione delle collisioni
- Simulazione della gravità
- Fisica realistica con rimbalzo e attrito

### PhysicsCatapultDemo
Gioco demo che dimostra le capacità del motore Artemis:
- **Gameplay**: Lancia proiettili da una catapulta per distruggere strutture
- **Fisica Realistica**: Tutti gli oggetti seguono le leggi della fisica
- **Materiali Diversi**: Blocchi in legno, pietra, vetro e metallo con proprietà fisiche uniche
- **Sistema di Danno**: I blocchi si danneggiano in base all'impatto
- **Tre Livelli**: Strutture diverse da abbattere

### AdvancedPhysicsDemo
Demo interattivo che mostra le funzionalità avanzate del motore Artemis:

**5 Scene Interattive:**

1. **Ragdoll Physics** - Sistema di ragdoll con articolazioni
   - RevoluteJoint e DistanceJoint per braccia e gambe
   - Fisica realistica del corpo umano
   - Spawn multipli di ragdoll

2. **Vehicle with Suspension** - Veicolo con sospensioni
   - SpringJoint per ammortizzatori realistici
   - Controllo accelerazione e salto
   - Terreno irregolare con ostacoli

3. **Chains & Bridges** - Catene e ponti sospesi
   - DistanceJoint per catene metalliche
   - RopeJoint per ponti flessibili
   - Fisica realistica della tensione

4. **Raycast Shooting** - Sistema di raycasting
   - Visualizzazione raggi in tempo reale
   - Layer mask filtering
   - Applicazione impulsi ai punti di impatto

5. **Trigger Zones & Events** - Zone trigger ed eventi
   - Collision events (OnEnter, OnStay, OnExit)
   - Trigger zones senza fisica
   - Monitoraggio eventi in tempo reale

## Controlli

- **W/S**: Regola l'angolo di lancio (su/giù)
- **A/D**: Regola la potenza del lancio (meno/più)
- **SPAZIO**: Lancia un proiettile
- **R**: Ricarica il livello corrente
- **1/2/3**: Carica il livello 1/2/3
- **Q**: Esci dal gioco

## Come Compilare ed Eseguire

```bash
# Compila l'intera solution
dotnet build Artemis.sln

# Esegui Physics Catapult Demo
dotnet run --project PhysicsCatapultDemo/PhysicsCatapultDemo.csproj

# Esegui Advanced Physics Demo
dotnet run --project AdvancedPhysicsDemo/AdvancedPhysicsDemo.csproj
```

## Caratteristiche del Motore Artemis

### Sistema Fisico Core
- Integrazione di Eulero per il movimento
- Rilevamento collisioni tra cerchi, rettangoli e forme miste
- Risoluzione delle collisioni con impulsi
- Attrito e rimbalzo configurabili
- Correzione della penetrazione
- Momento di inerzia per rotazioni realistiche
- Applicazione di forze, impulsi e torque

### Raycasting
- Sistema completo di raycast per cerchi e rettangoli
- Supporto per layer mask
- RaycastAll per rilevare tutti gli oggetti su un raggio
- Calcolo di punto di impatto, normale e distanza

### Sistema di Collisione Avanzato
- **Collision Layers**: Filtraggio bit-mask per controllare cosa collide con cosa
- **Collision Events**: Eventi OnCollisionEnter, OnCollisionStay, OnCollisionExit
- **Trigger Zones**: Sensori che rilevano collisioni senza fisica
- **Trigger Events**: OnTriggerEnter, OnTriggerStay, OnTriggerExit
- **Collision Listeners**: Interface per callback personalizzati

### Ottimizzazione Performance
- **Spatial Hash Grid**: Partizione spaziale per collisioni O(n) invece di O(n²)
- **QuadTree**: Alternativa al grid per distribuzioni non uniformi
- **Sleeping System**: Corpi inattivi vengono messi in "sleep" per risparmiare CPU
  - Wake up automatico su impulsi o collisioni
  - Soglie configurabili di velocità e tempo

### Joints e Constraints
- **DistanceJoint**: Mantiene distanza fissa tra due corpi
- **SpringJoint**: Molla con costante elastica e smorzamento
- **RevoluteJoint**: Rotazione attorno a un punto fisso (cerniera)
- **RopeJoint**: Vincolo di distanza massima
- Tutti i joint supportano stiffness configurabile

### Materiali
Ogni materiale ha proprietà uniche:

- **Legno**: Leggero, fragile, alta attrito
- **Pietra**: Pesante, resistente, molto attrito
- **Vetro**: Pesante ma fragile, basso attrito
- **Metallo**: Molto pesante, molto resistente

### Sistema di Rendering
Rendering ASCII/testuale per console che visualizza:
- Corpi fisici con caratteri diversi per materiale
- Indicatori di salute con cambio colore
- UI con statistiche di gioco
- Mirino della catapulta

## Architettura

Il progetto è diviso in due componenti principali:

1. **ArtemisEngine**: Libreria core del motore fisico
   - `Vector2.cs`: Matematica vettoriale 2D
   - `RigidBody.cs`: Corpi rigidi e forme con proprietà fisiche avanzate
   - `PhysicsWorld.cs`: Simulazione del mondo fisico con spatial partitioning
   - `Raycasting.cs`: Sistema di raycast per linee di vista e proiettili
   - `CollisionLayers.cs`: Sistema di layer per filtraggio collisioni
   - `CollisionEvents.cs`: Eventi e callback per collisioni
   - `SleepingSystem.cs`: Ottimizzazione per corpi inattivi
   - `Joints.cs`: Constraints e collegamenti tra corpi
   - `SpatialPartitioning.cs`: Spatial hash grid e quadtree

2. **PhysicsCatapultDemo**: Gioco dimostrativo
   - `Material.cs`: Definizione dei materiali
   - `GameObject.cs`: Oggetti di gioco con salute
   - `Catapult.cs`: Sistema di lancio
   - `LevelBuilder.cs`: Costruttore di livelli
   - `Renderer.cs`: Sistema di rendering ASCII
   - `Program.cs`: Game loop principale

3. **AdvancedPhysicsDemo**: Demo funzionalità avanzate
   - `RagdollBuilder.cs`: Creazione ragdoll con joints
   - `VehicleBuilder.cs`: Veicolo con sospensioni
   - `ChainBuilder.cs`: Catene e ponti sospesi
   - `AdvancedRenderer.cs`: Rendering con visualizzazione joints e raycast
   - `Program.cs`: 5 scene interattive

## Esempi di Utilizzo

### Creare un Corpo Rigido
```csharp
// Creare una sfera
var circleShape = new CircleShape(1.0f);
var ball = new RigidBody(new Vector2(10, 10), mass: 5.0f, circleShape);
ball.Restitution = 0.8f; // Rimbalzo
ball.Friction = 0.3f;

// Creare un rettangolo
var boxShape = new BoxShape(width: 2.0f, height: 1.0f);
var box = new RigidBody(new Vector2(5, 5), mass: 10.0f, boxShape);

// Aggiungere al mondo
world.AddBody(ball);
world.AddBody(box);
```

### Collision Layers
```csharp
// Impostare layer e mask
projectile.CollisionLayer = CollisionLayers.Projectile;
projectile.CollisionMask = CollisionLayers.Enemy | CollisionLayers.Ground;

enemy.CollisionLayer = CollisionLayers.Enemy;
enemy.CollisionMask = CollisionLayers.Projectile | CollisionLayers.Player;
```

### Raycasting
```csharp
var ray = new Ray(origin: new Vector2(0, 0),
                  direction: new Vector2(1, 0),
                  maxDistance: 100f);

var hit = world.Raycast(ray, layerMask: CollisionLayers.Enemy);
if (hit.Hit)
{
    Console.WriteLine($"Hit at {hit.Point}, distance: {hit.Distance}");
}
```

### Joints
```csharp
// Collegare due corpi con una molla
var spring = new SpringJoint(bodyA, bodyB,
                             restLength: 5.0f,
                             springConstant: 100f,
                             damping: 0.5f);
world.AddJoint(spring);

// Vincolo di distanza
var distance = new DistanceJoint(bodyA, bodyB, distance: 3.0f);
world.AddJoint(distance);
```

### Trigger Zones
```csharp
// Creare una zona trigger
var trigger = new RigidBody(new Vector2(20, 20), 0, new CircleShape(5f));
trigger.IsTrigger = true; // Non ha fisica, solo detection
trigger.IsStatic = true;

// Ascoltare eventi
world.OnTriggerEnter += (sender, e) =>
{
    Console.WriteLine($"Oggetto entrato nella trigger zone!");
};
```

### Collision Events
```csharp
world.OnCollisionEnter += (sender, e) =>
{
    Console.WriteLine($"Collisione tra {e.BodyA} e {e.BodyB}");
    // Applicare danno, suoni, effetti, ecc.
};
```

## Requisiti

- .NET 8.0 SDK o superiore
- Console con supporto per colori ANSI

## Licenza

Progetto dimostrativo per il motore fisico Artemis.
