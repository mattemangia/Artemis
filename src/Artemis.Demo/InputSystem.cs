using System;
using System.Collections.Generic;
using Artemis.Core;

namespace Artemis.Demo
{
    /// <summary>
    /// Abstract input system that supports mouse, touch, and keyboard input.
    /// Can be implemented for different platforms (Console, Unity, MonoGame, etc.).
    /// </summary>
    public interface IInputProvider
    {
        /// <summary>
        /// Gets the current pointer (mouse/touch) position in screen coordinates.
        /// </summary>
        Artemis.Demo.Vector2D PointerPosition { get; }

        /// <summary>
        /// Gets whether the primary pointer is pressed (left click / touch).
        /// </summary>
        bool IsPointerPressed { get; }

        /// <summary>
        /// Gets whether the secondary pointer is pressed (right click / two-finger tap).
        /// </summary>
        bool IsSecondaryPressed { get; }

        /// <summary>
        /// Gets all active touch points.
        /// </summary>
        IReadOnlyList<TouchPoint> TouchPoints { get; }

        /// <summary>
        /// Gets whether a key is currently pressed.
        /// </summary>
        bool IsKeyPressed(InputKey key);

        /// <summary>
        /// Gets whether a key was just pressed this frame.
        /// </summary>
        bool IsKeyJustPressed(InputKey key);

        /// <summary>
        /// Updates the input state. Call once per frame.
        /// </summary>
        void Update();
    }

    /// <summary>
    /// 2D vector for screen coordinates.
    /// </summary>
    public struct Vector2D
    {
        public double X;
        public double Y;

        public Vector2D(double x, double y)
        {
            X = x;
            Y = y;
        }

        public static Vector2D Zero => new(0, 0);

        public double Magnitude => Math.Sqrt(X * X + Y * Y);

        public static Vector2D operator +(Vector2D a, Vector2D b) => new(a.X + b.X, a.Y + b.Y);
        public static Vector2D operator -(Vector2D a, Vector2D b) => new(a.X - b.X, a.Y - b.Y);
        public static Vector2D operator *(Vector2D v, double s) => new(v.X * s, v.Y * s);

        /// <summary>
        /// Converts screen position to world position for a top-down view.
        /// </summary>
        public Vector3D ToWorld3D(
            double screenWidth,
            double screenHeight,
            double worldWidth,
            double worldHeight,
            double worldY = 0)
        {
            double worldX = (X / screenWidth - 0.5) * worldWidth;
            double worldZ = (Y / screenHeight - 0.5) * worldHeight;
            return new Vector3D(worldX, worldY, worldZ);
        }
    }

    /// <summary>
    /// Represents a touch point.
    /// </summary>
    public struct TouchPoint
    {
        public int Id;
        public Artemis.Demo.Vector2D Position;
        public TouchPhase Phase;
    }

    /// <summary>
    /// Touch phase (similar to Unity's TouchPhase).
    /// </summary>
    public enum TouchPhase
    {
        Began,
        Moved,
        Stationary,
        Ended,
        Canceled
    }

    /// <summary>
    /// Input keys enumeration.
    /// </summary>
    public enum InputKey
    {
        None,
        Left, Right, Up, Down,
        Space, Enter, Escape,
        A, B, C, D, E, F, G, H, I, J, K, L, M,
        N, O, P, Q, R, S, T, U, V, W, X, Y, Z,
        Num0, Num1, Num2, Num3, Num4, Num5, Num6, Num7, Num8, Num9
    }

    /// <summary>
    /// Console-based input provider.
    /// Uses arrow keys and WASD for movement, space for action.
    /// </summary>
    public class ConsoleInputProvider : IInputProvider
    {
        private Artemis.Demo.Vector2D _pointerPosition;
        private bool _isPointerPressed;
        private readonly HashSet<InputKey> _pressedKeys = new();
        private readonly HashSet<InputKey> _justPressedKeys = new();
        private readonly List<TouchPoint> _touchPoints = new();

        public Artemis.Demo.Vector2D PointerPosition => _pointerPosition;
        public bool IsPointerPressed => _isPointerPressed;
        public bool IsSecondaryPressed => false;
        public IReadOnlyList<TouchPoint> TouchPoints => _touchPoints;

        public double MoveSpeed { get; set; } = 0.5;
        public double ScreenWidth { get; set; } = 20;
        public double ScreenHeight { get; set; } = 20;

        public ConsoleInputProvider()
        {
            _pointerPosition = new Artemis.Demo.Vector2D(ScreenWidth / 2, ScreenHeight / 2);
        }

        public bool IsKeyPressed(InputKey key) => _pressedKeys.Contains(key);
        public bool IsKeyJustPressed(InputKey key) => _justPressedKeys.Contains(key);

        public void Update()
        {
            _justPressedKeys.Clear();
            _isPointerPressed = false;

            while (Console.KeyAvailable)
            {
                var keyInfo = Console.ReadKey(true);
                var key = MapKey(keyInfo.Key);

                if (key != InputKey.None)
                {
                    _pressedKeys.Add(key);
                    _justPressedKeys.Add(key);
                }

                // Update pointer position based on arrow keys
                switch (keyInfo.Key)
                {
                    case ConsoleKey.LeftArrow:
                    case ConsoleKey.A:
                        _pointerPosition.X = Math.Max(0, _pointerPosition.X - MoveSpeed);
                        break;
                    case ConsoleKey.RightArrow:
                    case ConsoleKey.D:
                        _pointerPosition.X = Math.Min(ScreenWidth, _pointerPosition.X + MoveSpeed);
                        break;
                    case ConsoleKey.UpArrow:
                    case ConsoleKey.W:
                        _pointerPosition.Y = Math.Max(0, _pointerPosition.Y - MoveSpeed);
                        break;
                    case ConsoleKey.DownArrow:
                    case ConsoleKey.S:
                        _pointerPosition.Y = Math.Min(ScreenHeight, _pointerPosition.Y + MoveSpeed);
                        break;
                    case ConsoleKey.Spacebar:
                    case ConsoleKey.Enter:
                        _isPointerPressed = true;
                        break;
                }
            }

            // Clear pressed keys after a frame (simplified)
            _pressedKeys.Clear();
        }

        private InputKey MapKey(ConsoleKey key) => key switch
        {
            ConsoleKey.LeftArrow => InputKey.Left,
            ConsoleKey.RightArrow => InputKey.Right,
            ConsoleKey.UpArrow => InputKey.Up,
            ConsoleKey.DownArrow => InputKey.Down,
            ConsoleKey.Spacebar => InputKey.Space,
            ConsoleKey.Enter => InputKey.Enter,
            ConsoleKey.Escape => InputKey.Escape,
            ConsoleKey.A => InputKey.A,
            ConsoleKey.B => InputKey.B,
            ConsoleKey.C => InputKey.C,
            ConsoleKey.D => InputKey.D,
            ConsoleKey.R => InputKey.R,
            ConsoleKey.Q => InputKey.Q,
            _ => InputKey.None
        };
    }

    /// <summary>
    /// Gesture recognizer for touch input.
    /// </summary>
    public class GestureRecognizer
    {
        private readonly List<TouchPoint> _activeTouches = new();
        private Artemis.Demo.Vector2D _lastTouchPosition;
        private DateTime _touchStartTime;
        private bool _isDragging;

        /// <summary>
        /// Minimum distance to be considered a drag.
        /// </summary>
        public double DragThreshold { get; set; } = 10;

        /// <summary>
        /// Maximum time for a tap gesture (milliseconds).
        /// </summary>
        public double TapTimeout { get; set; } = 300;

        /// <summary>
        /// Event raised when a tap is detected.
        /// </summary>
        public event Action<Artemis.Demo.Vector2D>? OnTap;

        /// <summary>
        /// Event raised when dragging starts.
        /// </summary>
        public event Action<Artemis.Demo.Vector2D>? OnDragStart;

        /// <summary>
        /// Event raised during dragging.
        /// </summary>
        public event Action<Artemis.Demo.Vector2D, Artemis.Demo.Vector2D>? OnDrag; // current, delta

        /// <summary>
        /// Event raised when dragging ends.
        /// </summary>
        public event Action<Artemis.Demo.Vector2D>? OnDragEnd;

        /// <summary>
        /// Event raised for pinch gesture (zoom).
        /// </summary>
        public event Action<double>? OnPinch; // scale factor

        /// <summary>
        /// Updates gesture recognition.
        /// </summary>
        public void Update(IReadOnlyList<TouchPoint> touches)
        {
            if (touches.Count == 0)
            {
                if (_isDragging)
                {
                    _isDragging = false;
                    OnDragEnd?.Invoke(_lastTouchPosition);
                }
                _activeTouches.Clear();
                return;
            }

            var touch = touches[0];

            switch (touch.Phase)
            {
                case TouchPhase.Began:
                    _lastTouchPosition = touch.Position;
                    _touchStartTime = DateTime.Now;
                    _isDragging = false;
                    break;

                case TouchPhase.Moved:
                    var delta = touch.Position - _lastTouchPosition;
                    if (!_isDragging && delta.Magnitude > DragThreshold)
                    {
                        _isDragging = true;
                        OnDragStart?.Invoke(touch.Position);
                    }

                    if (_isDragging)
                    {
                        OnDrag?.Invoke(touch.Position, delta);
                    }

                    _lastTouchPosition = touch.Position;
                    break;

                case TouchPhase.Ended:
                    if (!_isDragging)
                    {
                        var elapsed = (DateTime.Now - _touchStartTime).TotalMilliseconds;
                        if (elapsed < TapTimeout)
                        {
                            OnTap?.Invoke(touch.Position);
                        }
                    }
                    else
                    {
                        OnDragEnd?.Invoke(touch.Position);
                        _isDragging = false;
                    }
                    break;
            }

            // Handle pinch (two fingers)
            if (touches.Count >= 2)
            {
                var t1 = touches[0].Position;
                var t2 = touches[1].Position;
                double currentDistance = (t2 - t1).Magnitude;

                // Would need to track previous distance for proper pinch detection
                // This is a simplified version
            }
        }
    }

    /// <summary>
    /// Helper to convert screen coordinates to world physics coordinates.
    /// </summary>
    public class ScreenToWorldConverter
    {
        public double ScreenWidth { get; set; }
        public double ScreenHeight { get; set; }
        public AABB WorldBounds { get; set; }
        public double CameraHeight { get; set; } = 15;

        public ScreenToWorldConverter(double screenWidth, double screenHeight, AABB worldBounds)
        {
            ScreenWidth = screenWidth;
            ScreenHeight = screenHeight;
            WorldBounds = worldBounds;
        }

        /// <summary>
        /// Converts screen position to world position (top-down view).
        /// </summary>
        public Vector3D ScreenToWorld(Artemis.Demo.Vector2D screenPos, double worldY = 0)
        {
            double normalizedX = screenPos.X / ScreenWidth;
            double normalizedY = screenPos.Y / ScreenHeight;

            double worldX = WorldBounds.Min.X + normalizedX * (WorldBounds.Max.X - WorldBounds.Min.X);
            double worldZ = WorldBounds.Min.Z + normalizedY * (WorldBounds.Max.Z - WorldBounds.Min.Z);

            return new Vector3D(worldX, worldY, worldZ);
        }

        /// <summary>
        /// Converts world position to screen position.
        /// </summary>
        public Artemis.Demo.Vector2D WorldToScreen(Vector3D worldPos)
        {
            double normalizedX = (worldPos.X - WorldBounds.Min.X) / (WorldBounds.Max.X - WorldBounds.Min.X);
            double normalizedZ = (worldPos.Z - WorldBounds.Min.Z) / (WorldBounds.Max.Z - WorldBounds.Min.Z);

            return new Artemis.Demo.Vector2D(normalizedX * ScreenWidth, normalizedZ * ScreenHeight);
        }

        /// <summary>
        /// Casts a ray from screen position into the world.
        /// </summary>
        public (Vector3D origin, Vector3D direction) ScreenToRay(Artemis.Demo.Vector2D screenPos)
        {
            var worldPoint = ScreenToWorld(screenPos, 0);
            var origin = new Vector3D(worldPoint.X, CameraHeight, worldPoint.Z);
            var direction = Vector3D.Down;

            return (origin, direction);
        }
    }
}
