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

## Controlli

- **W/S**: Regola l'angolo di lancio (su/giù)
- **A/D**: Regola la potenza del lancio (meno/più)
- **SPAZIO**: Lancia un proiettile
- **R**: Ricarica il livello corrente
- **1/2/3**: Carica il livello 1/2/3
- **Q**: Esci dal gioco

## Come Compilare ed Eseguire

```bash
# Compila il progetto
dotnet build PhysicsCatapultDemo/PhysicsCatapultDemo.csproj

# Esegui il gioco
dotnet run --project PhysicsCatapultDemo/PhysicsCatapultDemo.csproj
```

## Caratteristiche del Motore Artemis

### Sistema Fisico
- Integrazione di Eulero per il movimento
- Rilevamento collisioni tra cerchi, rettangoli e forme miste
- Risoluzione delle collisioni con impulsi
- Attrito e rimbalzo configurabili
- Correzione della penetrazione

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
   - `RigidBody.cs`: Corpi rigidi e forme
   - `PhysicsWorld.cs`: Simulazione del mondo fisico

2. **PhysicsCatapultDemo**: Applicazione di gioco
   - `Material.cs`: Definizione dei materiali
   - `GameObject.cs`: Oggetti di gioco con salute
   - `Catapult.cs`: Sistema di lancio
   - `LevelBuilder.cs`: Costruttore di livelli
   - `Renderer.cs`: Sistema di rendering ASCII
   - `Program.cs`: Game loop principale

## Requisiti

- .NET 8.0 SDK o superiore
- Console con supporto per colori ANSI

## Licenza

Progetto dimostrativo per il motore fisico Artemis.
